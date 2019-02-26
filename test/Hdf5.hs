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

testData :: A.Vector Double
testData = A.fromList (A.Z A.:. 1000) [let x' = fromIntegral x in x' | x <- [0..]]

testDataC :: A.Vector ComplexDouble
testDataC = A.fromList (A.Z A.:. 1000) [let x' = fromIntegral x in x' :+ x' | x <- [0..]]

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
    remover filename
    createh5File filename
    createDatasetDouble filename datasetD testData
    createDatasetComplex filename datasetC testDataC
    testData' <- readDatasetDouble filename datasetD
    testDataC' <- readDatasetComplex filename datasetC
    let r1 = map (uncurry (==)) (zip (A.toList testData) (A.toList testData'))
    let r2 = map (uncurry (==)) (zip (A.toList testDataC) (A.toList testDataC'))
    remover filename
    return (and r1 && and r2)