{-# language TypeOperators       #-}
{-# language ViewPatterns        #-}
{-# language ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell #-}
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

import qualified Data.Array.Accelerate.LLVM.Native        as CPU
import qualified Data.Array.Accelerate.Interpreter        as I
import Data.Array.Accelerate.Debug                        as A

import qualified Prelude as P
import Prelude as P (fromIntegral, fromInteger, fromRational, String, return, (>>=), (>>), IO, Maybe(..), maybe)
import qualified Data.List as L (sort,(!!),sortBy)
import Data.Time.Clock

--aw_gridding :: String -> String -> String -> IO Image
aw_gridding :: Runners (String -> String -> String -> Maybe Int -> Maybe String -> IO (Scalar F))
aw_gridding run runN wfile afile datfile n outfile = do
    --setFlag dump_phases
    let theta    = 0.008
        lam      = 300000
    tim <- getCurrentTime
    P.putStrLn $ P.show tim
    vis <- readVis datfile 
    uvw <- readBaselines datfile
    (a1,a2,ts,f) <- readSource datfile
    let t = linearIndexArray ts 0
    akerns <- getAKernels afile theta t f
    wkerns <- getWKernels wfile theta
    let oargs = noOtherArgs
        akernels = use akerns
        wkernels = use (P.fst wkerns)
        wbins    = use (P.snd wkerns)
        args = noArgs
        len = constant . maybe (arraySize vis) P.id $ n
        len_ = lift (Z :. len) :: Exp DIM1
        src0 = zip4 (use a1) (use a2) (use ts) (fill len_ (constant f))
        uvwmat = use uvw
        u = slice uvwmat (constant (Z :. All :. (0 :: Int)))
        v = slice uvwmat (constant (Z :. All :. (1 :: Int)))
        w = slice uvwmat (constant (Z :. All :. (2 :: Int)))
        uvw0 = uvw_lambda f . take len $ zip3 u v w
        vis0 = take len $ use vis

        ones = fill len_ 1

        wt = doweight theta lam uvw0 ones
        (uvw1,vis1) = unzip $ mirror_uvw uvw0 vis0

        myuvw = uvw1
        mysrc = src0
        myvis = vis1
        
        {-
        uvgridf = $(CPU.runQ (aw_imaging noArgs noOtherArgs 0.008 300000))
        --uvgridf = runN (aw_imaging noArgs noOtherArgs 0.008 300000)
        uvgrid = uvgridf (CPU.run wkernels) (CPU.run wbins) (CPU.run akernels) (CPU.run myuvw) (CPU.run mysrc) (CPU.run $ zipWith (*) myvis wt)
        uvgrid1 = make_grid_hermitian (use uvgrid)
        -}
        uvgrid  = aw_imaging noArgs noOtherArgs 0.008 300000 wkernels wbins akernels myuvw mysrc (zipWith (*) myvis wt)
        uvgrid1 = make_grid_hermitian uvgrid

        img = run . map real . ifft $ uvgrid1
        max = run . maximum . flatten . use $ img
    P.putStrLn "Start imaging"
    case outfile of
        Nothing -> return ()
        Just fn -> do createh5File fn; createDatasetDouble fn "/img" img
    --return maxrun
    return $ max

    where
        readVis :: String -> IO (Vector Visibility)
        readVis file = do
            v <- readDatasetComplex file "/vis/vis" :: IO (Array DIM3 Visibility)
            let size = arraySize v
                sh   = arrayShape v
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


uvw_lambda :: Frequency -> Acc (Vector BaseLines) -> Acc (Vector BaseLines)
uvw_lambda f uvw = map mapper uvw
    where
        mapper :: Exp BaseLines -> Exp BaseLines
        mapper (unlift -> (u,v,w) :: (Exp F, Exp F, Exp F)) = lift (a * u,a *v, a * w)
        -- We divide by the speed of light and multiply with frequency
        a = constant $ f / 299792458.0