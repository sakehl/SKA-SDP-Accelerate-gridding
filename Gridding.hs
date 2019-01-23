{-# language FlexibleContexts    #-}
{-# language RebindableSyntax    #-}
{-# language ScopedTypeVariables #-}
{-# language TypeOperators       #-}
{-# language ViewPatterns        #-}


module Gridding where

import Data.Array.Accelerate                              as A hiding (fromInteger, fromRational, fromIntegral)
import qualified Data.Array.Accelerate                    as A (fromInteger, fromRational, fromIntegral)
import Data.Array.Accelerate.Data.Complex                 as A

import Data.Array.Accelerate.LLVM.Native                  as CPU
import Data.Array.Accelerate.Interpreter                  as I

import Data.Array.Accelerate.Math.DFT.Centre              as A
import Data.Array.Accelerate.Math.FFT                     as A

import qualified Prelude as P
import Prelude as P (fromIntegral, fromInteger, fromRational, String, return, (>>=), (>>), IO)
import Data.Char
import Text.Printf
import Data.List (intersperse, intercalate)

import Debug.Trace
import System.IO.Unsafe (unsafePerformIO)
import Control.Exception


-------------------------------
-- Gridding Accelerate code

-- # Simple imaging
simple_imaging :: Double                                 -- Field of view size
               -> Int                                    -- Grid size
               -> Vector (Double, Double, Double)        -- all the uvw baselines (coordinates) (lenght : n * (n-1))
               -> Vector (Int, Int, Double, Double)      -- NOT USED HERE (Antenna 1, Antenna 2, The time (in MJD UTC), Frequency (Hz)) (lenght : n * (n-1))
               -> Vector (Complex Double)                -- visibility  (length n * (n-1))
               -> Acc (Matrix (Complex Double))
simple_imaging theta lam uvw src vis =
    let lamf = fromIntegral lam
        n = P.round(theta * lamf)
        ne = constant n
        p = map (`div3` constant lamf) (use uvw)
        czero = constant $ 0 :+ 0
        guv = fill (index2 ne ne) czero
    in grid guv p (use vis)    

grid :: Acc (Matrix (Complex Double))                    -- Destination grid N x N of complex numbers
     -> Acc (Vector (Double, Double, Double))            -- The uvw baselines, but between -.5 and .5
     -> Acc (Vector (Complex Double))                    -- The visiblities
     -> Acc (Matrix (Complex Double))
grid a p v = permute (+) a (\i -> xy ! i) v
    where
        Z :. n :. _ = unlift (shape a) :: Z :. Exp Int :. Exp Int
        halfn = n `div` 2
        nf = A.fromIntegral n :: Exp Double
        xy = map gridxy p :: Acc (Vector DIM2)

        -- Note the order, python specified as a[y, x] += v
        -- We do the same
        gridxy :: Exp (Double, Double, Double) -> Exp DIM2
        gridxy (unlift -> (x, y, _ :: Exp Double)) = index2 (toGridCell y) (toGridCell x)

        toGridCell :: Exp Double -> Exp Int
        toGridCell f = halfn + floor (0.5 + nf * f)

-- # Normal imaging
frac_coord :: Exp Int                                   -- The height/width
           -> Exp Int                                   -- Oversampling factor
           -> Acc (Vector Double)                       -- a uvw baseline, but between -.5 and .5
           -> Acc (Vector (Int, Int))
frac_coord n qpx p = 
    let halfn = n `div` 2
        halfnf = A.fromIntegral halfn
        nf = A.fromIntegral n
        qpxf = A.fromIntegral qpx
        qpxfrac = 0.5 / qpxf

        x   = map (\a -> halfnf + a * nf) p
        flx = map (\a -> floor (a + qpxfrac) ) x
        fracx = zipWith (\a b -> round ((a - (A.fromIntegral b) )* qpxf) ) x flx
    in zip flx fracx

frac_coords :: Exp (Int, Int)                           -- (Heigt, width) of the grid
            -> Exp Int                                  -- Oversampling factor
            -> Acc (Vector (Double, Double, Double))    -- The uvw baselines, but between -.5 and .5
            -> Acc (Vector (Int, Int, Int, Int))
-- NOTE we should give it in height, width order.
frac_coords (unlift -> (h, w) ) qpx p =
    let (u, v, _) = unzip3 p
        (x, xf)   = unzip (frac_coord w qpx u)
        (y, yf)   = unzip (frac_coord h qpx v)
    in zip4 x xf y yf

convgrid :: Array DIM4 (Complex Double)                  -- The oversampled convolution kernel
         -> Acc (Matrix (Complex Double))                -- Destination grid N x N of complex numbers
         -> Acc (Vector (Double, Double, Double))        -- The uvw baselines, but between -.5 and .5
         -> Acc (Vector (Complex Double))                -- The visiblities
         -> Acc (Matrix (Complex Double))
convgrid gcf a p v =
    let Z :. qpx :. _ :. gh :. gw = (arrayShape gcf) :: Z :. Int :. Int :. Int :. Int
        Z :. height :. width = unlift (shape a) :: Z :. Exp Int :. Exp Int
        coords = frac_coords (lift (height, width)) (constant qpx) p
        gpugcf = use gcf
        halfgh = constant $ gh `div` 2
        halfgw = constant $ gw `div` 2

        getComplex :: Int -> Int -> Exp (Int, Int, Int, Int) -> Exp (Int, Int, Complex Double)
        getComplex i j (unlift -> (x, xf, y, yf)::(Exp Int,Exp Int,Exp Int,Exp Int)) =
            lift ( x
                 , y
                 , gpugcf ! lift (Z :. yf :. xf :. constant i :. constant j) )
        fullrepli :: Int -> Int -> Acc (Matrix (Complex Double)) -> Acc (Matrix (Complex Double))
        fullrepli i j origin = 
            let temp = map (getComplex i j) coords
                (x, y, source) = unzip3 temp
                indexer id = lift (Z :. y ! id - halfgh + constant i :. x ! id - halfgw + constant j)
                newV = zipWith (*) source v
            in permute (+) origin indexer newV
        
    in myfor gh (\i -> myfor gw (fullrepli i)) a

myfor :: Int -> (Int -> a -> a) -> a -> a
myfor n f x | n P.== 0  = x
            | P.otherwise =  myfor (n-1) f (f (n-1) x)


type ImagingFunction = Double                                  -- Field of view size
                     -> Int                                    -- Grid size
                     -> Vector (Double, Double, Double)        -- all the uvw baselines (coordinates) (lenght : n * (n-1))
                     -> Vector (Int, Int, Double, Double)      -- NOT USED HERE (Antenna 1, Antenna 2, The time (in MJD UTC), Frequency (Hz)) (lenght : n * (n-1))
                     -> Vector (Complex Double)                -- visibility  (length n * (n-1))
                     -> Acc (Matrix (Complex Double))

do_imaging :: Double                                 -- Field of view size
           -> Int                                    -- Grid size
           -> Vector (Double, Double, Double)        -- all the uvw baselines (coordinates) (lenght : n * (n-1))
           -> Vector (Int, Int, Double, Double)      -- NOT USED HERE (Antenna 1, Antenna 2, The time (in MJD UTC), Frequency (Hz)) (lenght : n * (n-1))
           -> Vector (Complex Double)                -- visibility  (length n * (n-1))
           -> ImagingFunction
           -> Acc (Matrix (Complex Double))
do_imaging theta lam uvw src vis imgfn = undefined

mirror_uvw :: Acc (Vector (Double, Double, Double)) -> Acc (Vector (Complex Double)) -> Acc (Vector ((Double, Double, Double), Complex Double))
mirror_uvw uvw vis = let
        (u,v,w) = unzip3 uvw
        uvwvis = zip4 u v w vis
        mirrorv o@(unlift -> (u,v,w,vis) :: (Exp Double, Exp Double, Exp Double, Exp (Complex Double))) =
            if v < 0 
                then lift ((u, -v, w), conjugate vis)
                else lift ((u,v,w),vis)
    in map mirrorv uvwvis


doweight :: Double -> Int -> Acc (Vector (Double, Double, Double)) -> Acc (Vector (Complex Double)) -> Acc (Vector (Complex Double))
doweight theta lam p v = let
        n = P.round(theta * lamf)
        lamf = fromIntegral lam
        ne = constant n
        lamfe = constant lamf
        coords = frac_coords (lift (ne, ne)) 1 (map (`div3` lamfe) p)

    in undefined
----------------------
-- Fourier transformations
ishift2D :: Elt a => Acc (Matrix a) -> Acc (Matrix a)
ishift2D = shift2D

fft :: Acc (Matrix (Complex Double)) -> Acc (Matrix (Complex Double))
fft = shift2D . fft2D Forward . ishift2D

ifft :: Acc (Matrix (Complex Double)) -> Acc (Matrix (Complex Double))
ifft = shift2D . fft2D Inverse . ishift2D

------------------------
-- Helper functions
div3 :: Exp (Double, Double, Double) -> Exp Double -> Exp (Double, Double, Double)
div3 (unlift -> (a,b,c)) x = lift (a / x, b / x, c / x)

----------------------
-- Testing stuff
testData :: ([Complex Double], [(Double, Double, Double)])
testData = unsafePerformIO testData0

testSimple :: Int -> Int -> IO ()
testSimple i j = do
    testData <- testData0
    let n   = (P.length . P.fst) testData
        m   = (P.length . P.snd) testData
        size = assert (n P.== m) $ Z :. n
        vis = fromList size $ P.fst testData
        uvw = fromList size $ P.snd testData
        src      = undefined
        theta    = 2*0.05
        lam      = 18000
        result = simple_imaging theta lam uvw src vis
        runresult = CPU.run result
        fst = indexArray runresult (Z :. i :. j)
        maxi = P.maximum (toList runresult)
    --P.writeFile "result.csv" (makeFile runresult)
    P.putStrLn (P.show maxi)

instance (P.Ord a) => P.Ord (Complex a)
    where
        (x1 :+ y1) <=  (x2 :+ y2) | x1 P.== x2 = y1 P.<= y2
                                  | P.otherwise = x1 P.<= x2

makeFile :: Matrix (Complex Double) -> String
makeFile mat = let
    (Z :. n :. m) = arrayShape  mat
    ls =  P.map showC $ toList mat
    lss = toRows n [] ls :: [[String]]
    rows = P.map (intercalate ",") lss :: [String]

    toRows _ res [] = P.reverse res
    toRows n res ys = let (x, xs) = P.splitAt n ys
                      in toRows n (x : res) xs

    showC (x :+ y) | y P.>= 0 = '(' : printf "%e" x P.++ '+' : printf "%e" y P.++ "j)"
                   | P.otherwise = '(' : printf "%e" x P.++ printf "%e" y P.++ "j)"
    in P.unlines rows


-------------------------
-- Test input reading
testData0 :: IO ([Complex Double], [(Double, Double, Double)])
testData0 = do
    visData <- P.readFile "vis.csv"
    uvwData <- P.readFile "uvw.csv"
    let vislines = P.lines visData
        visdat   = P.map parseComplex vislines
        uvwlines = P.lines uvwData
        uvwdat   = P.map parseTuple3 uvwlines
    return (visdat, uvwdat)

parseComplex :: String -> Complex Double
parseComplex s  =
    case P.head s of
        ' ' -> parseComplex (P.tail s)
        '(' -> parseComplex (P.tail s)
        _   -> let (x, rest) = readDouble s
                   (y, _)    = readDouble rest
               in P.read x :+ P.read y
    where
        readDouble s@('-':_) = P.splitAt 25 s
        readDouble ('+':s)   = P.splitAt 24 s
        readDouble s         = P.splitAt 24 s

parseTuple3 :: P.Read a => String -> (a,a,a)
parseTuple3 s = 
    let valid c = c P./= ',' P.&& P.not (isSpace c)
        (x, _ : rest0) = P.span valid s
        (y, _ : rest1) = P.span valid rest0
        (z, _)         = P.span valid rest1
    in (P.read x, P.read y, P.read z)

testing :: Acc (Vector Int)
testing = use $ fromList (Z :. 5) [0..]

reptesting :: Acc (Array DIM3 Int)
reptesting = replicate (constant (Z :. All :. (3::Int) :. (3:: Int))) $ testing

testing2 :: Exp Int -> Exp Int -> Exp Int -> (Exp Int)
testing2 a b c = reptesting ! index3 a b c

runExp :: Elt e => Exp e -> e
runExp e = indexArray (CPU.run (unit e)) Z

testing0 :: Vector (Int)
testing0 = fromList (Z :. 5) [0..]

testing3 :: Int
testing3 = indexArray testing0 (Z :. 4)

ffttest0 :: Int -> (Vector Double, Vector Double, Vector Double, [(Double, Int)])
ffttest0 n = let
    halfn = fromIntegral $ (n - 1) `div` 2
    minstart = case P.even n of
        True -> -halfn - 1
        False -> -halfn
    lst = [0..halfn] P.++ [minstart..(-1)]
    freqs = fromList (Z :. n) lst

    res  =  shift1D' (use freqs)
    res2 = ishift1D' (use freqs)
    res3 = shiftTest lst
    in (freqs, CPU.run res,CPU.run res2, res3)


ishift1D' :: Elt e => Acc (Vector e) -> Acc (Vector e)
ishift1D' arr
  = A.backpermute (A.shape arr) p arr
  where
    p ix
      = let Z:.x = unlift ix :: Z :. Exp Int
        in index1 (x A.< mw2 ? (x + mw, x - mw2))
    Z:.w    = unlift (A.shape arr)
    mw      = w `div` 2
    mw2     = (w + 1) `div` 2

shift1D' :: Elt e => Acc (Vector e) -> Acc (Vector e)
shift1D' arr
  = A.backpermute (A.shape arr) p arr
  where
    p ix
      = let Z:.x = unlift ix :: Z :. Exp Int
        in index1 (x A.< mw2 ? (x + mw, x - mw2))
    Z:.w    = unlift (A.shape arr)
    mw      = w `div` 2
    mw2     = (w + 1) `div` 2

shiftTest :: [Double] -> [(Double, Int)]
shiftTest xs = let
    temp = addindex 0 xs
    addindex n [] = []
    addindex n (x:xs) = (x, n) : addindex (n+1) xs

    w = P.length xs
    mw      = w `P.div` 2
    mw2     = (w + 1) `P.div` 2
    indexer (v,x) = case x P.< mw2 of
        True  -> (v, x + mw)
        False -> (v, x -mw2)
    in P.map indexer temp

--runAcc = CPU.run 
{-
def convgrid(
    -- The oversampled convolution kernel (most of the time it is shaped like 2 x 2 x 31 x 31)
    gcf, :: 
    --Destination grid N x N of complex numbers
    a,  :: Matrix (Complex Double)
    -- The uvw baselines, but between -.5 and .5
    p,  :: [(Double, Double, Double)]
    -- The visiblities
    v :: [Complex Double]
    ):
    """Grid after convolving with gcf

    Takes into account fractional `uv` coordinate values where the GCF
    is oversampled

    :param a: Grid to add to
    :param p: UVW positions
    :param v: Visibility values
    :param gcf: Oversampled convolution kernel
    """

    Qpx, _, gh, gw = gcf.shape
    coords = frac_coords(a.shape, Qpx, p)
    for v, x,xf, y,yf in zip(v, *coords):
        a[y-gh//2 : y+(gh+1)//2,
          x-gw//2 : x+(gw+1)//2] += gcf[yf,xf] * v
-}

-------------------------------
-- Gridding Python code
{-
def simple_imaging(
    --Field of view size
    theta, :: Double
    --Grid size
    lam,  :: Int
    -- all the uvw baselines (coordinates) (lenght : n * (n-1))                          
    uvw,  :: [(Double, Double, Double)]
    -- (Antenna 1, Antenna 2, The time (in MJD UTC), Frequency (Hz)) (lenght : n * (n-1))                     
    src,  :: [(Int, Int, Double, Double)] 
    -- The visibility 
    vis :: [Complex Double]
    ):
    """Trivial function for imaging
    Does no convolution but simply puts the visibilities into a grid
    cell i.e. boxcar gridding
    """
    N = int(round(theta * lam))
    assert N > 1
    //The uv plane to grid to N x N matrix of complex numbers
    guv = numpy.zeros([N, N], dtype=complex)
    grid(guv, uvw / lam, vis)
return guv


def grid(
    --Destination grid N x N of complex numbers
    a, :: Matrix (Complex Double)
    -- The uvw baselines, but between -.5 and .5
    p, :: [(Double, Double, Double)]
    -- The visiblities
    v :: [Complex Double]
    ):
    """Grid visibilities (v) at positions (p) into (a) without convolution
    :param a:   The uv plane to grid to (updated in-place!)
    :param p:   The coordinates to grid to (in fraction [-.5,.5[ of grid)
    :param v:   Visibilities to grid
    """
    assert numpy.max(p) < 0.5

    N = a.shape[0]
    xy = N//2 + numpy.floor(0.5 + N * p[:,0:2]).astype(int)
    for (x, y), v in zip(xy, v):
        a[y, x] += v

def frac_coord(
    -- The height/width
    N, :: Int
    -- Oversampling factor
    Qpx, :: Int
    -- a uvw baseline, but between -.5 and .5
    p, :: [Double]
    ):
    """
    Compute whole and fractional parts of coordinates, rounded to
    Qpx-th fraction of pixel size

    The fractional values are rounded to nearest 1/Qpx pixel value. At
    fractional values greater than (Qpx-0.5)/Qpx coordinates are
    roundeded to next integer index.

    :param N: Number of pixels in total
    :param Qpx: Fractional values to round to
    :param p: Coordinate in range [-.5,.5[
    """
    assert (p >= -0.5).all() and (p < 0.5).all()
    x = N//2 + p * N
    flx = numpy.floor(x + 0.5 / Qpx)
    fracx = numpy.around((x - flx) * Qpx)
    return flx.astype(int), fracx.astype(int)

def frac_coords(shape, Qpx, p):
    """Compute grid coordinates and fractional values for convolutional
    gridding

    :param shape: (height,width) grid shape
    :param Qpx: Oversampling factor
    :param p: array of (x,y) coordinates in range [-.5,.5[
    """
    h, w = shape # NB order (height,width) to match numpy!
    x, xf = frac_coord(w, Qpx, p[:,0])
    y, yf = frac_coord(h, Qpx, p[:,1])
    return x,xf, y,yf

def convgrid(
    -- The oversampled convolution kernel (most of the time it is shaped like 2 x 2 x 31 x 31)
    gcf, :: 
    --Destination grid N x N of complex numbers
    a,  :: Matrix (Complex Double)
    -- The uvw baselines, but between -.5 and .5
    p,  :: [(Double, Double, Double)]
    -- The visiblities
    v :: [Complex Double]
    ):
    """Grid after convolving with gcf

    Takes into account fractional `uv` coordinate values where the GCF
    is oversampled

    :param a: Grid to add to
    :param p: UVW positions
    :param v: Visibility values
    :param gcf: Oversampled convolution kernel
    """

    Qpx, _, gh, gw = gcf.shape
    coords = frac_coords(a.shape, Qpx, p)
    for v, x,xf, y,yf in zip(v, *coords):
        a[y-gh//2 : y+(gh+1)//2,
          x-gw//2 : x+(gw+1)//2] += gcf[yf,xf] * v
-}