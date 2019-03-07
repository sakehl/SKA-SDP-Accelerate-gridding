{-# language TypeOperators       #-}
module ImageDataset where

import Types
import Gridding
import Hdf5
import Text.Printf

import Data.Array.Accelerate                              as A hiding (fromInteger, fromRational, fromIntegral)
import qualified Data.Array.Accelerate                    as A (fromInteger, fromRational, fromIntegral)
import Data.Array.Accelerate.Data.Complex                 as A
import Data.Array.Accelerate.Math.DFT.Centre              as A
import Data.Array.Accelerate.Math.FFT                     as A

import Data.Array.Accelerate.LLVM.Native                  as CPU
import Data.Array.Accelerate.Interpreter                  as I

import qualified Prelude as P
import Prelude as P (fromIntegral, fromInteger, fromRational, String, return, (>>=), (>>), IO, Maybe(..), maybe)
import qualified Data.List as L (sort,(!!) )

aw_gridding :: String -> String -> String -> IO Image
aw_gridding wfile afile datfile = do
    let theta    = 0.08
        lam      = 180
    vis <- readVis datfile 
    uvw <- readBaselines datfile
    (a1,a2,ts,f) <- readSource datfile
    let t = linearIndexArray ts 0
    akerns <- getAKernels afile theta t f
    wkerns <- getWKernels wfile
    let oargs = noOtherArgs{akernels = Just akerns, wkernels = Just wkerns}
        args = noArgs
        (res, _, _) =CPU.run $ do_imaging theta lam uvw a1 a2 ts f vis (aw_imaging args oargs)
        --testing = CPU.run $ aw_imaging args oargs theta lam (use uvw) (use src) (use vis)
    return res 

    where
        readVis :: String -> IO (Vector Visibility)
        readVis file = readDatasetComplex file "vis/vis"

        readBaselines :: String -> IO (Matrix BaseLine)
        readBaselines file = readDatasetDouble file "/vis/uvw"

        readSource :: String -> IO (Vector Antenna, Vector Antenna, Vector Time, Frequency)
        readSource file = do
            a1 <- readDatasetInt file "/vis/antenna1" :: IO (Vector Antenna)
            a2 <- readDatasetInt file "/vis/antenna2" :: IO (Vector Antenna)
            t  <- readDatasetDouble file "/vis/time" :: IO (Vector Time)
            f  <- readDatasetDouble file "/vis/frequency" :: IO (Vector Frequency)
            let f0 = linearIndexArray f 0
            return (a1, a2, t, f0)

        getAKernels :: String -> F -> Time -> Frequency -> IO (Array DIM3 Visibility)
        getAKernels file theta t f = do 
            ants <- listGroupMembers file (printf "/akern/%f" theta)
            
            
            
            return (fromList (Z :. 3 :. 4 :. 4) [ sin (5 * x + 1) :+ cos x| x <- [0..] ])

        getWKernels :: String -> IO (Array DIM5 Visibility, Vector BaseLine)
        getWKernels _ = return (fromList (Z :. 2 :. 2 :. 2 :. 4 :. 4) [ sin (5 * x + 1) :+ cos x| x <- [0..] ], fromList (Z :. 2) [0..])
        
getAKernels :: String -> F -> Time -> Frequency -> IO (Array DIM3 Visibility)
getAKernels file theta t f = do
    P.putStrLn "Loading A kernels"
    --Get all the antennas
    ants <- listGroupMembers file (printf "/akern/%f" theta)
    let ants0 = L.sort . P.map P.read $ ants :: [Int]
        a0 = P.head ants0
    --Get all the times
    ts <- listGroupMembers file (printf "/akern/%f/%i" theta a0)
    let ts0 = L.sort . P.map P.read $ ts :: [Time]
        closestt = findClosestL ts0 t
    P.putStrLn $ printf "Closest time found in a-kernels is %f (input time: %f)" closestt t
    --Get all the frequencies
    fs <- listGroupMembers file (printf "/akern/%f/%i/%.1f" theta a0 closestt)
    let fs0 = L.sort . P.map P.read $ fs :: [Frequency]
        closestf = findClosestL fs0 f
    P.putStrLn $ printf "Closest frequency found in a-kernels is %f (input frequency: %f)" closestf f
    
    -- For now only get the first one
    akern0 <- readDatasetComplex file (printf "/akern/%f/%i/%.1f/%.0f/kern" theta a0 closestt closestf) :: IO (Matrix Visibility)
    allkerns <- P.mapM (\a -> readDatasetComplex file (printf "/akern/%f/%i/%.1f/%.0f/kern" theta a closestt closestf)) ants0 :: IO [Matrix Visibility]
    let Z :. y :. x = arrayShape akern0 :: DIM2
        newsh = Z :. 1 :. y :. x :: DIM3
        akern0' = arrayReshape newsh akern0
    return akern0'

test = getAKernels "SKA1_Low_akern.h5" 0.08 0 0

findClosestL :: (P.Num a, P.Ord a) => [a] -> a -> a
findClosestL ws w =
    let ifte g i e | g = i
                   | P.otherwise = e
        cmp (min, max) = (max P.- min) `P.div` 2 P.>= 1
        f (min, max) = 
            let id = (max P.+ min) `P.div` 2
            in ifte (w P.> ws L.!! id)
                (id, max)
                (min, id)
        max = P.length ws - 1
        minmax = (0, max)
        
        (r1, r2) = mywhile cmp f $ minmax
        v1 = ws L.!! r1
        v2 = ws L.!! r2

    in ifte (P.abs (w P.- v1) P.< abs (w P.- v2)) v1 v2

mywhile :: (a -> Bool) -> (a -> a) -> a -> a
mywhile cmp f init | P.not (cmp init) = init
                 | P.otherwise = mywhile cmp f (f P.$! init)

concatMatrices :: Elt a => [Matrix a] -> Array DIM3 a
concatMatrices ms = undefined