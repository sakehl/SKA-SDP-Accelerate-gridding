{-# LANGUAGE RebindableSyntax #-}
{-# LANGUAGE TypeApplications #-}

import Hdf5

import qualified Data.Array.Accelerate             as A
import Data.Array.Accelerate.Debug       as A
import Data.Array.Accelerate.Data.Bits   as A
import Data.Array.Accelerate.Interpreter as I
import Data.Array.Accelerate.LLVM.Native as CPU
import Data.Array.Accelerate.Data.Complex as A

import qualified Data.Vector.Storable as SV 

import Prelude hiding (catch)
import System.Directory
import Control.Exception
import System.IO.Error hiding (catch)
import Foreign.C.Types

testData :: A.Matrix Double
testData = A.fromList (A.Z A.:. 1000 A.:. 2) [let x' = fromIntegral x in x' | x <- [0..]]

testDataC :: A.Vector ComplexDouble
testDataC = A.fromList (A.Z A.:. 1000) [let x' = fromIntegral x in x' :+ x' | x <- [0..]]

test3 :: A.Array A.DIM4 ComplexDouble
test3 = A.fromList (A.Z A.:. 2 A.:. 3 A.:. 1000 A.:. 2) [let x' = fromIntegral x in x' :+ x' | x <- [0..]]

remover :: FilePath -> IO ()
remover fileName = removeFile fileName `catch` handle
    where handle e
            | isDoesNotExistError e = return ()
            | otherwise = throwIO e

testIO :: IO Bool
testIO = do
    let filename = "test.h5"
    let datasetD = "\testD"
    let datasetC = "\testC"
    let dataset3 = "\test3"
    remover filename
    createh5File filename
    createDatasetDouble filename datasetD testData
    createDatasetComplex filename datasetC testDataC
    createDatasetComplex filename dataset3 test3
    testData' <- readDatasetDouble filename datasetD
    testDataC' <- readDatasetComplex filename datasetC
    test3' <- readDatasetComplex filename dataset3
    --let comparer = CPU.runN (\xs ys -> A.zipWith (A.==) xs ys)
    let r1 = comparer testData testData' `A.indexArray` A.Z
    let r2 = comparer testDataC testDataC' `A.indexArray` A.Z
    let r3 = comparer test3 test3' `A.indexArray` A.Z
    remover filename
    return (r1 && r2 && r3)
        where
            comparer :: (A.Elt e, A.Shape sh, A.Eq e) => A.Array sh e -> A.Array sh e -> A.Scalar Bool
            comparer = CPU.runN (\xs ys -> A.foldAll (A.&&) (A.lift True) $ A.zipWith (A.==) xs ys)