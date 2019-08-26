{-# language FlexibleContexts    #-}
{-# language ScopedTypeVariables #-}
{-# language TypeOperators       #-}
{-# language ViewPatterns        #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE RankNTypes          #-}

module Main where

import Data.Array.Accelerate                              as A hiding (fromInteger, fromRational, fromIntegral)
import qualified Data.Array.Accelerate                    as A (fromInteger, fromRational, fromIntegral)
import Data.Array.Accelerate.Data.Complex                 as A
import Data.Array.Accelerate.Math.DFT.Centre              as A
import Data.Array.Accelerate.Math.FFT                     as A
-- import Data.Array.Accelerate.Math.FFT.Type                as A

import Data.Array.Accelerate.LLVM.Native                  as CPU
import Data.Array.Accelerate.Interpreter                  as I
import Data.Array.Accelerate.Debug                        as A hiding (Mode)

import Data.Array.Accelerate.Array.Sugar                  as S hiding (shape)
import Data.Array.Accelerate.LLVM.Native.Foreign
import Data.Array.Accelerate.Type                         as A
import Math.FFT                                           as FFT
import Math.FFT.Base                                      ( FFTWReal, Sign(..), Flag, measure, destroyInput )
import Data.Array.CArray                                            ( CArray )

import Gridding
import Types
import Hdf5
import ImageDataset
import FFT
import Utility

import qualified Prelude as P
import Prelude as P (fromIntegral, fromInteger, fromRational, String, return, (>>=), (>>), IO, Maybe(..), maybe)
import Data.Char (isSpace)
import Text.Printf
import Data.List (intercalate)
import Data.Time.Clock

import Debug.Trace
import System.IO.Unsafe (unsafePerformIO)
import Control.Exception (assert)

import Control.Lens as L (_1, _2, (^.) )

main :: IO ()
main = do 
    let n_         = Just 2
        dir        = "data"
        wkern      = dir P.++ "/" P.++ "SKA1_Low_wkern2.h5"
        akern      = dir P.++ "/" P.++ "SKA1_Low_akern3.h5"
        visdata    = dir P.++ "/" P.++ "SKA1_Low_quick.h5"
    
    --f <- Main.aw_gridding CPU.run wkern akern visdata n_
    --P.print f
    --setflags
    P.print . CPU.run $ testLoop1
    P.print . CPU.run $ testLoop2
    P.print . CPU.run $ testLoop3
    P.print . CPU.run $ testLoop4
    return ()

xs' :: Vector Int
xs' = fromList (Z :. 12) [0..]

xs2' :: Vector Int
xs2' = fromList (Z :. 9) [10..]

xs :: Acc (Vector Int)
xs = use xs'

ys :: Matrix Int
ys = fromList (Z :. 2 :. 10) [0..]

ysSeq1 :: Seq [Vector Int]
ysSeq1 = toSeqInner . use $ ys

ysSeq2 :: Seq [Vector Int]
ysSeq2 = toSeqOuter . use $ ys

ysSeq3 :: Seq [Matrix Int]
ysSeq3 = subarrays (lift $ (Z :. 1 :. 10 :: Z :. Int :. Int)) $ ys

ysSeq4 :: Seq [Vector Int]
ysSeq4 = streamIn [xs', xs2', xs', xs']

check :: Shape sh => Acc (Array sh Int) -> Acc (Scalar Bool) 
check = map (<100) . sum

checkF :: Shape sh => Acc (Array sh Int) -> Acc (Scalar Bool) 
checkF _ = unit . constant $ False

check2 :: Shape sh => Acc (Array sh Int, Scalar Int) -> Acc (Scalar Bool) 
check2 = map (<10) . asnd

f :: Shape sh => Acc (Array sh Int) -> Acc (Array sh Int)
f = map (+1)

f2 :: Shape sh => Acc (Array sh Int, Scalar Int) -> Acc (Array sh Int, Scalar Int)
f2 (unlift -> (a,n) :: (Acc (Array sh Int), Acc (Scalar Int))) = lift (map (+1) a, map (+1) n)

fSeq :: Shape sh => Seq [Array sh Int] -> Seq [Array sh Int]
fSeq = mapSeq f

loop :: Shape sh => Acc (Array sh Int) -> Acc (Array sh Int)
loop = awhile check f

loopSeq :: Shape sh => Seq [Array sh Int] -> Seq [Array sh Int]
loopSeq = mapSeq loop

testLoop1 = collect . tabulate . loopSeq $ ysSeq1
testLoop2 = collect . tabulate . loopSeq $ ysSeq2
testLoop3 = collect . tabulate . loopSeq $ ysSeq3
testLoop4 = collect . tabulate . loopSeq $ ysSeq4


flags = [dump_phases]

setflags :: IO ()
setflags = do
    setFlag dump_cc
    setFlag dump_sched
    setFlags Main.flags

unsetflags :: IO ()
unsetflags = clearFlags Main.flags

----------------------
-- I/O things (for testing)

-- Extensive img data set things
aw_gridding :: Runners0 -> String -> String -> String -> Maybe Int -> IO (Scalar F)
aw_gridding run wfile afile datfile n = do
    --setFlag dump_phases
    let theta    = 0.008
        lam      = 300000
    tim <- getCurrentTime
    P.putStrLn $ P.show tim
    vis <- readVis datfile 
    uvw <- readBaselines datfile
    (a1,a2,ts,f) <- readSource datfile
    let t = linearIndexArray ts 0
        f' = use $ fromList Z [f]
    akerns <- getAKernels afile theta t f False
    wkerns <- getWKernels wfile theta False
    let oargs = noOtherArgs{akernels = Just akerns, wkernels = Just wkerns}
        args = noArgs
        akernels = use akerns
        wkernels = use (P.fst wkerns)
        wbins    = use (P.snd wkerns)
        --(res, _, _) =run $ do_imaging theta lam uvw a1 a2 ts f vis (aw_imaging args oargs)
        len = constant . maybe (arraySize . arrayShape $ vis) P.id $ n
        len_ = lift (Z :. len) :: Exp DIM1
        --len = arrayShape vis
        src0 = zip4 (use a1) (use a2) (use ts) (fill len_ (constant f))
        uvwmat = use uvw
        u = slice uvwmat (constant (Z :. All :. (0 :: Int)))
        v = slice uvwmat (constant (Z :. All :. (1 :: Int)))
        w = slice uvwmat (constant (Z :. All :. (2 :: Int)))
        uvw0 = uvw_lambda f' . take len $ zip3 u v w
        vis0 = take len $ use vis

        ones = fill len_ 1
        wt = doweight theta lam uvw0 ones
        (uvw1,vis1) = unzip $ mirror_uvw uvw0 vis0

        myuvw = uvw1
        mysrc = src0
        myvis = vis1
        uvgrid = aw_imaging args oargs theta lam wkernels wbins akernels myuvw mysrc myvis
        maxuvgridrun = run . maximum . flatten . map real $ uvgrid

        uvgrid1 = make_grid_hermitian uvgrid

        img = map real . shift2D . myfft2D Inverse . ishift2D $ uvgrid1
        
        max = (maximum . flatten) img
        (imgrun, maxrun) = run $ (lift (img,max) :: Acc (Matrix F, Scalar F))

        akern1 = slice (use akerns) (constant (Z :. (0 :: Int) :. All :. All))
        akern2 = slice (use akerns) (constant (Z :. (1 :: Int) :. All :. All))
        wkern  = slice (use (P.fst wkerns)) (constant (Z :. (0 :: Int) :. (0 :: Int):. (0 :: Int):. All :. All))
        aakern = convolve2d akern1 akern2
        awkern = convolve2d wkern aakern
 
        lamf = fromIntegral lam
        nn = constant . P.round $ theta * lamf

        p = map (`div3` constant lamf) myuvw
        guv = fill (index2 nn nn) 0 :: Acc (Matrix Visibility)

        (_,_,ww) = unzip3 myuvw
        closestw = map (findClosest (wbins)) ww
        (a11, a22, _, _) = unzip4 mysrc
        index = zipWith3 (\x y z -> lift (x, A.fromIntegral y, A.fromIntegral z)) closestw a11 a22 :: Acc (Vector (Int, Int, Int))

        Z :. _ :. qpx :. _ :. gh :. gw = unlift . shape . use . P.fst $ wkerns :: Z :. Exp Int :. Exp Int :. Exp Int :. Exp Int :. Exp Int
        Z :. height :. width = unlift (shape guv) :: Z :. Exp Int :. Exp Int
        halfgh = gh `div` 2
        halfgw = gw `div` 2
        coords = frac_coords (lift (height, width)) qpx p
        (cx, cxf, cy, cyf) = unzip4 coords
        dat = zip5 cx cxf cy cyf myvis

        (unlift -> (x, xf, y, yf, vis_) :: (Exp Int, Exp Int, Exp Int, Exp Int, Exp Visibility) ) = dat !! 0 

        (unlift -> (wbin, a1index, a2index)  :: (Exp Int, Exp Int, Exp Int)) = index !! 0 
        
        a1_ = slice (use akerns) (lift (Z :. a1index :. All :. All))
        a2_ = slice (use akerns) (lift (Z :. a2index :. All :. All))
        w_  = slice (use (P.fst wkerns)) (lift (Z :. wbin :. All :. All :. All :. All))
        
        aakern2 = convolve2d a1_ a2_
        awkern2 = map conjugate $ aw_kernel_fn2 yf xf w_ a1_ a2_
        

    --P.putStrLn (P.show $ run $ minimum . map real . flatten $ uvgrid)
    P.putStrLn "Start imaging"
    --createh5File "result.h5"
    --createDatasetDouble "result.h5" "/img" imgrun
    {-
    createh5File "convolveTest.h5"
    createDatasetComplex "convolveTest.h5" "/a0" (run $ akern1)
    createDatasetComplex "convolveTest.h5" "/a1" (run $ akern2)
    createDatasetComplex "convolveTest.h5" "/w" (run $ wkern)
    createDatasetComplex "convolveTest.h5" "/aa" (run $ aakern)
    createDatasetComplex "convolveTest.h5" "/aw" (run $ awkern)
    createDatasetComplex "convolveTest.h5" "/uvgrid" (run $ uvgrid)
    createDatasetComplex "convolveTest.h5" "/awkern" (run $ awkern2)
    createDatasetComplex "convolveTest.h5" "/herm" (run $ uvgrid1)
    createDatasetDouble "convolveTest.h5" "/img" imgrun
    P.putStrLn $ printf "myuvw: %s\n mysrc: %s \n mvis: %s" (P.show . CPU.run . unit . shape $ u) (P.show . CPU.run . unit . shape $ src0) (P.show . CPU.run $ myvis)
    --P.mapM P.putStrLn . P.map P.show . toList . CPU.run $ uvw1
    --P.putStrLn . P.show . CPU.run $ unit $ lift (wbins !! (closestw !! 0), a11 !! 0, a22 !! 0)
    --P.putStrLn . P.show . CPU.run $ coords
    --P.putStrLn . P.show . CPU.run $ myvis
    -}
    return maxuvgridrun
    --return $ fromList Z [1.0]

    where
        readVis :: String -> IO (Vector Visibility)
        readVis file = do
            v <- readDatasetComplex file "/vis/vis" :: IO (Array DIM3 Visibility)
            let size = arraySize sh
                sh   = arrayShape v
                newv = arrayReshape (Z :. size :: DIM1) v
            P.putStrLn . P.show $ sh 
            return newv

        readBaselines :: String -> IO (Matrix BaseLine)
        readBaselines file = do
            r <- readDatasetDouble file "/vis/uvw"
            return r

        readSource :: String -> IO (Vector Antenna, Vector Antenna, Vector Time, Frequency)
        readSource file = do
            a1 <- readDatasetInt64 file "/vis/antenna1" :: IO (Vector Antenna)
            a2 <- readDatasetInt64 file "/vis/antenna2" :: IO (Vector Antenna)
            t  <- readDatasetDouble file "/vis/time" :: IO (Vector Time)
            f  <- readDatasetDouble file "/vis/frequency" :: IO (Vector Frequency)
            let f0 = linearIndexArray f 0
            return (a1, a2, t, f0)

        linearIndexArray :: Vector a -> Int -> a
        linearIndexArray a n = toList a P.!! n



-- myconvolve2d :: Acc (Matrix Visibility) -> Acc (Matrix Visibility) -> Acc (Matrix Visibility)
-- myconvolve2d a1 a2 = 
--     let
--         (Z :. n :. _) = (unlift . shape) a1 :: Z :. Exp Int :. Exp Int
--         m_ = 2*n - 1
--         -- We need m to be a power of 2
--         pw = A.ceiling (logBase 2 (A.fromIntegral m_) :: Exp F) :: Exp Int
--         m = 32 --2 ^ pw
--         m2 = A.fromIntegral $ m * m
        

--         fft mode = if True then myfft2DAv mode else fft2D mode
        
--         a1fft = fft Inverse . ishift2D $ pad_mid a1 m
--         a2fft = fft Inverse . ishift2D $ pad_mid a2 m
--         a1a2 = zipWith (*) a1fft a2fft
--         convolved = shift2D . fft Forward $ a1a2
--         mid = extract_mid n convolved
--     in map (*m2) mid


-- fft2DAvoid :: forall e. (Numeric e) 
--                     => Mode 
--                     -> ForeignAcc (Array DIM2 (Complex e) -> Array DIM2 (Complex e))
-- fft2DAvoid mode = ForeignAcc (nameOf mode (undefined::DIM2))
--     $ case numericR::NumericR e of
--               NumericRfloat32 -> liftIO . liftAtoC go
--               NumericRfloat64 -> liftIO . liftAtoC go
--     where
--         go :: FFTWReal r => CArray (Int,Int) (Complex r) -> CArray (Int,Int) (Complex r)
--         go = FFT.dftG (signOf mode) fftflags [0,1]

-- myfft2DAv :: Numeric e => Mode -> Acc (Array DIM2 (Complex e)) -> Acc (Array DIM2 (Complex e))
-- myfft2DAv mode = foreignAcc (fft2DAvoid mode) $ A.map (\_ -> -89978978.4e0100) {-Bollocks implementation, for checking-}

-- fftflags :: Flag
-- fftflags = estimate .|. destroyInput