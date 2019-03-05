module ImageDataset where

import Types
import Gridding

import Data.Array.Accelerate                              as A hiding (fromInteger, fromRational, fromIntegral)
import qualified Data.Array.Accelerate                    as A (fromInteger, fromRational, fromIntegral)
import Data.Array.Accelerate.Data.Complex                 as A
import Data.Array.Accelerate.Math.DFT.Centre              as A
import Data.Array.Accelerate.Math.FFT                     as A

import Data.Array.Accelerate.LLVM.Native                  as CPU
import Data.Array.Accelerate.Interpreter                  as I

import qualified Prelude as P
import Prelude as P (fromIntegral, fromInteger, fromRational, String, return, (>>=), (>>), IO, Maybe(..), maybe)

aw_gridding :: String -> String -> String -> IO (Image)
aw_gridding wfile afile datfile = do
    vis <- readVis datfile 
    uvw <- readBaselines datfile
    src <- readSource datfile
    let f = undefined
        t = undefined
    akerns <- getAKernels afile t f
    wkerns <- getWKernels wfile
    let oargs = noOtherArgs{akernels = Just akerns, wkernels = Just wkerns}
        args  = noArgs
        theta    = 0.08
        lam      = 180
        (res, _, _) =CPU.run $ do_imaging theta lam uvw src vis (aw_imaging args oargs)
    return res

    where
        readVis :: String -> IO (Vector Visibility)
        readVis _ = return (fromList (Z :. 2) [x :+ 0.5 | x <- [0..] ])

        readBaselines :: String -> IO (Vector BaseLines)
        readBaselines _ = return (fromList (Z :. 2) [(x*2, x*3, x) | x <- [0..] ])

        readSource :: String -> IO (Vector (Antenna, Antenna, Time, Frequency))
        readSource _ = return (fromList (Z :. 2) [(x, x + 1, 0 , 0) | x <- [0..] ])

        readFreq :: String -> IO Frequency
        readFreq _ = return 5

        getAKernels :: String -> Time -> Frequency -> IO (Array DIM3 Visibility)
        getAKernels _ _ _ = return (fromList (Z :. 3 :. 4 :. 4) [ sin (5 * x + 1) :+ cos x| x <- [0..] ])

        getWKernels :: String -> IO (Array DIM5 Visibility, Vector BaseLine)
        getWKernels _ = return (fromList (Z :. 2 :. 2 :. 2 :. 4 :. 4) [ sin (5 * x + 1) :+ cos x| x <- [0..] ], fromList (Z :. 2) [0..])
        