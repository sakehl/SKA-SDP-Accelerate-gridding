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
import Data.Array.Accelerate.Debug                        as A

import qualified Prelude as P
import Prelude as P (fromIntegral, fromInteger, fromRational, String, return, (>>=), (>>), IO, Maybe(..), maybe)
import qualified Data.List as L (sort,(!!),sortBy)

--aw_gridding :: String -> String -> String -> IO Image
aw_gridding :: String -> String -> String -> IO FourierSpace
aw_gridding wfile afile datfile = do
    setFlag dump_phases
    let theta    = 0.08
        lam      = 180
    vis <- readVis datfile 
    uvw <- readBaselines datfile
    (a1,a2,ts,f) <- readSource datfile
    let t = linearIndexArray ts 0
    akerns <- getAKernels afile theta t f
    wkerns <- getWKernels wfile theta
    let oargs = noOtherArgs{akernels = Just akerns, wkernels = Just wkerns}
        args = noArgs
        (res, _, _) =CPU.run $ do_imaging theta lam uvw a1 a2 ts f vis (aw_imaging args oargs)
        len = arrayShape vis
        src0 = zip4 (use a1) (use a2) (use ts) (fill (lift len) (constant f))
        uvwmat = use uvw
        u = slice uvwmat (constant (Z :. All :. (0 :: Int)))
        v = slice uvwmat (constant (Z :. All :. (1 :: Int)))
        w = slice uvwmat (constant (Z :. All :. (2 :: Int)))
        uvw0 = zip3 u v w
        myuvw = (uvw0)
        mysrc = (src0)
        myvis = (use vis)
        testing = CPU.run $ aw_imaging args oargs theta lam myuvw mysrc myvis
    P.putStrLn "Start imaging"
    --P.putStrLn $ printf "myuvw: %s\n mysrc: %s \n mvis: %s" (P.show . CPU.run . unit . shape $ u) (P.show . CPU.run . unit . shape $ src0) (P.show . CPU.run $ myvis)
    return testing

    where
        readVis :: String -> IO (Vector Visibility)
        readVis file = do
            v <- readDatasetComplex file "/vis/vis" :: IO (Array DIM3 Visibility)
            let size = arraySize v
                newv = arrayReshape (Z :. size :: DIM1) v
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
        
getAKernels :: String -> F -> Time -> Frequency -> IO (Array DIM3 Visibility)
getAKernels file theta t f = do
    P.putStrLn "Loading A kernels"
    --Get all the antennas
    ants <- listGroupMembers file (printf "/akern/%f" theta)
    let antsSorted = convertAndSort ants :: [(Int, String)]
        a0 = P.snd . P.head $ antsSorted
    --Get all the times
    ts <- listGroupMembers file (printf "/akern/%f/%s" theta a0)
    let tsSorted = convertAndSort ts :: [(Time, String)]
        ts0      = P.map P.fst tsSorted
        idt      = P.snd $ findClosestList ts0 t
        closestt =  P.snd $ tsSorted L.!! idt
    P.putStrLn $ printf "Closest time found in a-kernels is %s (input time: %f)" closestt t
    --Get all the frequencies
    fs <- listGroupMembers file (printf "/akern/%f/%s/%s" theta a0 closestt)
    let fsSorted = convertAndSort fs :: [(Frequency, String)]
        fs0      = P.map P.fst tsSorted
        idf      = P.snd $ findClosestList fs0 f
        closestf =  P.snd $ fsSorted L.!! idf
    P.putStrLn $ printf "Closest frequency found in a-kernels is %s (input frequency: %f)" closestf f
    -- Get all the kernels from the antenna's
    let alldatasets = P.map  (\a -> printf "/akern/%f/%s/%s/%s/kern" theta (P.snd a) closestt closestf) antsSorted
    akerns <- readDatasetsComplex file alldatasets :: IO (Array DIM3 Visibility)
    P.putStrLn "A kernels loaded"
    return akerns

    
getWKernels :: String -> F -> IO (Array DIM5 Visibility, Vector BaseLine)
getWKernels file theta = do
    P.putStrLn "Loading W kernels"
    --Get all the wbins
    wbins <- listGroupMembers file (printf "/wkern/%f" theta)
    let wbinsSorted = convertAndSort wbins :: [(BaseLine, String)]
        alldatasets = P.map  (\w -> printf "/wkern/%f/%s/kern" theta (P.snd w)) wbinsSorted
        --The wbins as vector
        wbinsv = fromList (Z :. P.length wbinsSorted :: DIM1) (P.map (P.fst) wbinsSorted)
    -- Get all the kernels from the wbins
    wkerns <- readDatasetsComplex file alldatasets :: IO (Array DIM5 Visibility)
    P.putStrLn "W kernels loaded"
    return (wkerns, wbinsv)




-- Some helpers
findClosestList :: (P.Num a, P.Ord a) => [a] -> a -> (a, Int)
findClosestList ws w =
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

    in ifte (P.abs (w P.- v1) P.< abs (w P.- v2)) (v1, r1) (v2, r2)

mywhile :: (a -> Bool) -> (a -> a) -> a -> a
mywhile cmp f init | P.not (cmp init) = init
                   | P.otherwise = mywhile cmp f (f P.$! init)

convertAndSort :: (P.Read a, P.Ord a) => [String] -> [(a, String)]
convertAndSort xs = let
    converted = P.map (\x -> (P.read x, x) ) xs
    sorter (a,_) (b,_) = P.compare a b
    in L.sortBy sorter converted