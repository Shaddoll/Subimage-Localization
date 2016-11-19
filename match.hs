module Main where
import qualified Vision.Image as I
import Vision.Image.Storage.DevIL
import Vision.Primitive
import qualified Data.Vector.Storable as V
import Data.Array.CArray
import qualified Data.Array.IArray as Arr
import qualified Foreign.Storable as ST
import Math.FFT
import Math.FFT.Base
import Data.Complex
import Data.Word
import Control.Parallel.Strategies
import Control.Exception.Base
import Control.DeepSeq
import System.Environment (getArgs)

getHeight :: DIM2 -> Int
getHeight ((Z :. h) :. w) = h

getWidth :: DIM2 -> Int
getWidth ((Z :. h) :. w) = w

meanFilter' :: (Int, Int) -> I.Grey -> I.Manifest (Double)
meanFilter' (bw, bh) grey = I.Manifest (I.manifestSize grey) (V.fromList aft)
  where
    w = getWidth (I.manifestSize grey)
    h = getHeight (I.manifestSize grey)
    raw = vectorToMatrix (I.manifestVector grey) (h - 1, w - 1)
    aft = boxAverageFilter raw (w, h) (bw, bh)

matrixToVector :: CArray (Int, Int) Double -> V.Vector Double
matrixToVector m = V.fromList $ elems m

getGreyPixel :: I.GreyPixel -> Double
getGreyPixel (I.GreyPixel x) = fromIntegral x

subGrey :: I.GreyPixel -> Double -> Double
subGrey (I.GreyPixel x) y = fromIntegral x - y

threadDouble :: Double -> Double
threadDouble x
  | y >= 0.001 = y
  | otherwise = 0.001
    where
      y = sqrt $ abs x

threadPixel :: Double -> Double -> Double
threadPixel x y
  | z > 127 = 127
  | z < -127 = -127
  | otherwise = z
    where
      z = x * 32 / y

listToMatrix :: ST.Storable a => [a] -> (Int, Int) -> CArray (Int, Int) a
listToMatrix !l !(x, y) = Arr.listArray ((0, 0), (x, y)) l

vectorToMatrix ::  V.Vector I.GreyPixel -> (Int, Int) -> CArray (Int, Int) Double
vectorToMatrix v (x, y) = listToMatrix (V.toList (V.map getGreyPixel v)) (x, y)

vectorToMatrix' ::  V.Vector Double -> (Int, Int) -> CArray (Int, Int) Double
vectorToMatrix' v (x, y) = listToMatrix (V.toList v) (x, y)

extendMatrixUnsafe :: (ST.Storable a, Num a) => CArray (Int, Int) a -> (Int, Int) -> CArray (Int, Int) a
extendMatrixUnsafe matrix (x , y) = array ((0, 0), (x, y)) ((assocs matrix) ++ [((i, j), 0) | i <- [upperX + 1 .. x], j <- [upperY + 1 .. y]])
    where
      upperX = fst $ snd (bounds matrix)
      upperY = snd $ snd (bounds matrix)

conjugateComplexMatrix :: (FFTWReal r) => CArray (Int, Int) (Complex r) -> CArray (Int, Int) (Complex r)
conjugateComplexMatrix a = amap conjugate a

crossMultiplyComplexMatrix :: (FFTWReal r) => CArray (Int, Int) (Complex r) -> CArray (Int, Int) (Complex r) -> CArray (Int, Int) (Complex r)
crossMultiplyComplexMatrix a b = array ((0, 0), (upperX, upperY)) [((i, j), (a ! (i, j) * (b ! (i, j)))) | i <- [0 .. upperX], j <- [0 .. upperY]]
    where
      upperX = fst $ snd (bounds a)
      upperY = snd $ snd (bounds a)

fft2 :: (FFTWReal r) => CArray (Int, Int) r -> CArray (Int, Int) (Complex r)
fft2 = dftRCN [0, 1]

ifft2 :: (FFTWReal r) => CArray (Int, Int) (Complex r) -> CArray (Int, Int) r
ifft2 = dftCRN [0, 1]

rgbToGrey :: I.RGB -> I.Grey
rgbToGrey = I.convert

loadBMPRGB :: FilePath -> IO (Either StorageError I.RGB)
loadBMPRGB = load BMP

maxMatrixIndex :: (FFTWReal r) => CArray (Int, Int) r -> (Int, Int)
maxMatrixIndex matrix = maxArrayIndexP matrix (0, 0) (snd (bounds matrix)) (-1, -1) (-1 / 0)
maxArrayIndexP :: (FFTWReal r) => CArray (Int, Int) r -> (Int, Int) -> (Int, Int) -> (Int, Int) -> r -> (Int, Int)
maxArrayIndexP matrix (x, y) (upperX, upperY) result num
  | y > upperY = result
  | x < upperX - 1 = if (current > num) then maxArrayIndexP matrix (x + 1, y) (upperX, upperY) (x, y) current else maxArrayIndexP matrix (x + 1, y) (upperX, upperY) result num
  | otherwise = if (current > num) then maxArrayIndexP matrix (0, y + 1) (upperX, upperY) (x, y) current else maxArrayIndexP matrix (0, y + 1) (upperX, upperY) result num
    where
      current = matrix ! (x, y)

meanFilter :: (Int, Int) -> I.Grey -> I.Manifest (Double)
meanFilter (w, h) = I.mean ((Z :. h) :. w)

round' :: Double -> Word16
round' = round

meanFilterVector :: (Int, Int) -> (Int, Int) -> V.Vector Double -> V.Vector Double
meanFilterVector (bw, bh) (w, h) v = I.manifestVector(I.mean ((Z :. bh) :. bw) (I.Manifest ((Z :. h) :. w) (V.map round' v)))

exchange :: (a, b) -> (b, a)
exchange (x, y) = (y, x)

main :: IO ()
main = do
  [big, part, mode] <- getArgs
  Right img_partial <- loadBMPRGB part
  Right img <- loadBMPRGB big
  let (Z :. p_height) :. p_width = I.manifestSize img_partial
  let (Z :. height) :. width = I.manifestSize img
  (normalPartial, normal) <- evaluate (runEval (parallelProcess img img_partial (width, height) (p_width, p_height) (p_width `div` 2, p_height `div` 2)))
  let corr = ifft2 (crossMultiplyComplexMatrix normal normalPartial)
  let r = maxMatrixIndex corr
  print $ exchange r

processImagePartial :: I.RGB -> (Int, Int) -> (Int, Int) -> (Int, Int) -> CArray (Int, Int) (Complex Double)
processImagePartial rgb (w, h) (bw, bh) (nw, nh) = result
  where
    grey = rgbToGrey rgb
    grey_v = I.manifestVector grey
    meanfilter = meanFilter' (bw, bh) grey
    mean_v = I.manifestVector meanfilter
    diff = V.zipWith subGrey grey_v mean_v
    m = extendMatrixUnsafe (vectorToMatrix' diff (h - 1, w - 1)) (nh - 1, nw - 1)
    result = conjugateComplexMatrix (fft2 m)

processImage :: I.RGB -> (Int, Int) -> (Int, Int) -> CArray (Int, Int) (Complex Double)
processImage rgb (w, h) (bw, bh) = result
  where
    grey = rgbToGrey rgb
    grey_v = I.manifestVector grey
    meanfilter = meanFilter' (bw, bh) grey
    mean_v = I.manifestVector meanfilter
    diff = V.zipWith subGrey grey_v mean_v
    m = vectorToMatrix' diff (h - 1, w - 1)
    result = fft2 m

parallelProcess :: I.RGB -> I.RGB -> (Int, Int) -> (Int, Int) -> (Int, Int) -> Eval (CArray (Int, Int) (Complex Double), CArray (Int, Int) (Complex Double))
parallelProcess rgb rgb_p (w, h) (pw, ph) (bw, bh) = do
  b <- rpar (processImage rgb (w, h) (bw, bh))
  a <- rseq (processImagePartial rgb_p (pw, ph) (bw, bh) (w, h))
  b <- rseq b
  return (a, b)

boxAverageFilter :: CArray (Int, Int) Double -> (Int, Int) -> (Int, Int) -> [Double]
boxAverageFilter matrix img_size@(img_width, img_height) (box_width, box_height) = boxVerticalAverageFilter (listToMatrix (boxHorizontalAverageFilter matrix img_size box_width) (img_height-1,img_width-1)) img_size box_height

boxHorizontalAverageFilter :: CArray (Int, Int) Double -> (Int, Int) -> Int -> [Double]
boxHorizontalAverageFilter !matrix !img_size !box_width = boxHorizontalAverageFilterP matrix img_size box_width (box_width - box_width `div` 2) (0, 0) (0, 0) []
boxHorizontalAverageFilterP :: CArray (Int, Int) Double -> (Int, Int) -> Int -> Int -> (Int, Int) -> (Double, Int) -> [Double] -> [Double]
boxHorizontalAverageFilterP matrix img_size@(img_width, img_height) box_width init_count !(x,y) !(sum,count) !result
  | y >= img_height = force result
  | x == 0 = boxHorizontalAverageFilterP matrix img_size box_width init_count (x+1,y) (start_sum,init_count) ((start_sum / fromIntegral init_count):result)
  | x < 1 + box_width `div` 2 = boxHorizontalAverageFilterP matrix img_size box_width init_count (x+1,y) (left_sum,count+1) ((left_sum / fromIntegral(count+1)):result)
  | x < img_width + 1 - box_width + box_width `div` 2 = boxHorizontalAverageFilterP matrix img_size box_width init_count (x+1,y) (middle_sum,count) ((middle_sum / fromIntegral count):result)
  | x < img_width - 1 = boxHorizontalAverageFilterP matrix img_size box_width init_count (x+1,y) (right_sum, count-1) ((right_sum / fromIntegral (count-1)):result)
  | otherwise = boxHorizontalAverageFilterP matrix img_size box_width init_count (0,y+1) (right_sum, count-1) ((right_sum / fromIntegral (count-1)):result)
    where
      start_sum = foldl (+) 0 (map (matrix ! ) [(y, j) | j <- [0 .. init_count - 1]])
      left_sum = sum + matrix ! (y, x + (box_width - box_width `div` 2 - 1))
      middle_sum = sum + (matrix ! (y, x + (box_width - box_width `div` 2 - 1)) - matrix ! (y, x - (box_width `div` 2 - 1)))
      right_sum = sum - (matrix ! (y, x - (box_width `div` 2 - 1)))

boxVerticalAverageFilter :: CArray (Int, Int) Double -> (Int, Int) -> Int -> [Double]
boxVerticalAverageFilter matrix img_size box_height = boxHorizontalAverageFilterP matrix img_size box_height (box_height - box_height `div` 2) (0, 0) (0, 0) []
boxVerticalAverageFilterP :: CArray (Int, Int) Double -> (Int, Int) -> Int -> Int -> (Int, Int) -> (Double, Int) -> [Double] -> [Double]
boxVerticalAverageFilterP matrix img_size@(img_width, img_height) box_height init_count !(x,y) !(sum,count) !result
  | x >= img_width = force result
  | y == 0 = boxVerticalAverageFilterP matrix img_size box_height init_count (x,y+1) (start_sum,init_count) ((start_sum / fromIntegral init_count):result)
  | y < 1 + box_height `div` 2 = boxVerticalAverageFilterP matrix img_size box_height init_count (x,y+1) (bottom_sum,count+1) ((bottom_sum / fromIntegral(count+1)):result)
  | y < img_height + 1 - box_height + box_height `div` 2 = boxVerticalAverageFilterP matrix img_size box_height init_count (x,y+1) (middle_sum,count) ((middle_sum / fromIntegral count):result)
  | y < img_height - 1 = boxVerticalAverageFilterP matrix img_size box_height init_count (x,y+1) (top_sum, count-1) ((top_sum / fromIntegral (count-1)):result)
  | otherwise = boxVerticalAverageFilterP matrix img_size box_height init_count (x+1,0) (top_sum, count-1) ((top_sum / fromIntegral (count-1)):result)
    where
      start_sum = foldl (+) 0 (map (matrix ! ) [(i, x) | i <- [0 .. init_count - 1]])
      bottom_sum = sum + (matrix ! (y + (box_height - box_height `div` 2 - 1), x))
      middle_sum = sum + (matrix ! (y + (box_height - box_height `div` 2 - 1), x) - matrix ! (y, x - (box_height `div` 2 - 1)))
      top_sum = sum - (matrix ! (y - (box_height `div` 2 - 1), x))
