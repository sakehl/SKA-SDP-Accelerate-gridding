{-# LANGUAGE ForeignFunctionInterface #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE TypeFamilies #-}

module Hdf5 where

import Foreign
import Foreign.C.Types
import Foreign.C.String
import Foreign.Ptr
import System.IO.Unsafe (unsafePerformIO)
import Text.Printf
import Control.Exception (assert)
import System.Directory
import Control.Monad

import qualified Data.Vector.Storable as SV
import qualified Data.Array.Accelerate.IO.Data.Vector.Storable as A
import qualified Data.Array.Accelerate                         as A
import qualified Data.Array.Accelerate.IO.Foreign.ForeignPtr   as A
import qualified Data.Array.Accelerate.Array.Sugar             as A hiding(shape)
import Data.Array.Accelerate.Data.Complex

type ComplexCDouble = Complex CDouble
type ComplexDouble = Complex Double

---- Actual functionallity
foreign import ccall unsafe "createh5File"
    c_createh5File :: CString -> IO ()

foreign import ccall unsafe "getRankDataset"
    c_getRankDataset :: CString -> CString -> IO CInt

foreign import ccall unsafe "getDimsDataset"
    c_getDimsDataset :: CString -> CString -> CInt -> Ptr CInt -> IO ()

foreign import ccall unsafe "readDatasetInt"
    c_readDatasetInt :: CString -> CString -> Ptr CLong -> IO ()

foreign import ccall unsafe "readDatasetLLong"
    c_readDatasetInt64 :: CString -> CString -> Ptr CLong -> IO ()

foreign import ccall unsafe "readDatasetDouble"
    c_readDatasetDouble :: CString -> CString -> Ptr CDouble -> IO ()

foreign import ccall unsafe "readDatasetComplex"
    c_readDatasetComplex :: CString -> CString -> Ptr ComplexCDouble -> IO ()

foreign import ccall unsafe "readDatasetsComplex"
    c_readDatasetsComplex :: CString -> Ptr CString -> Ptr ComplexCDouble -> IO ()

foreign import ccall unsafe "readDatasetsDouble"
    c_readDatasetsDouble :: CString -> Ptr CString -> Ptr CDouble -> IO ()

foreign import ccall unsafe "createDatasetInt"
    c_createDatasetInt64 :: CString -> CString -> CInt -> Ptr CInt -> Ptr CLong -> IO ()

foreign import ccall unsafe "createDatasetDouble"
    c_createDatasetDouble :: CString -> CString -> CInt -> Ptr CInt -> Ptr CDouble -> IO ()

foreign import ccall unsafe "createDatasetComplex"
    c_createDatasetComplex :: CString -> CString -> CInt -> Ptr CInt -> Ptr ComplexCDouble -> IO ()

foreign import ccall unsafe "listGroupMembers"
    c_listGroupMembers :: CString -> CString -> IO (Ptr CString)

createh5File :: String -> IO ()
createh5File name' = do 
    name <- newCString name'
    c_createh5File name

createDatasetDoubleV :: String -> String -> CInt -> SV.Vector CInt -> SV.Vector CDouble -> IO ()
createDatasetDoubleV = createDatasetV c_createDatasetDouble
createDatasetComplexV :: String -> String -> CInt -> SV.Vector CInt -> SV.Vector ComplexCDouble -> IO ()
createDatasetComplexV = createDatasetV c_createDatasetComplex

createDatasetV :: Storable a => (CString -> CString -> CInt -> Ptr CInt -> Ptr a -> IO ()) -> 
    String -> String -> CInt -> SV.Vector CInt -> SV.Vector a -> IO ()
createDatasetV c_createDataset name' dataset' rank dims dat = do
    name <- newCString name'
    dataset <- newCString dataset'
    SV.unsafeWith dims $ \dimsp -> SV.unsafeWith dat $ \datap -> c_createDataset name dataset rank dimsp datap

getDimsDataset :: String -> String -> CInt -> IO (SV.Vector CInt)
getDimsDataset name' dataset' rank = do
    name <- newCString name'
    dataset <- newCString dataset'
    -- Allocate enough memory for the array
    dims <- mallocArray (fromIntegral rank)
    -- Get the dimensions via the FFI binding
    c_getDimsDataset name dataset rank dims
    -- Make it to a foreign poitner, this will make sure it gets freed eventually and we can get a Vector out of it
    foreign_dims <- newForeignPtr finalizerFree dims
    -- Give us back a vector :)
    return $ SV.unsafeFromForeignPtr0 foreign_dims (fromIntegral rank)

getRankDataset :: String -> String -> IO CInt
getRankDataset name' dataset' = do
    name <- newCString name'
    dataset <- newCString dataset'
    c_getRankDataset name dataset

listGroupMembers :: String -> String -> IO [String]
listGroupMembers name' group' = do
    name <- newCString name'
    group <- newCString group'
    listp <- c_listGroupMembers name group
    cstringp <- peekArray0 nullPtr listp
    mapM peekCString cstringp

readDataset :: forall e a b sh.(Storable b, A.EltRepr b ~ e, A.Elt b, A.ForeignPtrs e ~ ForeignPtr a, A.Shape sh) =>
             (CString -> CString -> Ptr b -> IO ())
            -> String -> String -> IO (A.Array sh b)
readDataset reader name' dataset' = do
    exists <- doesFileExist name'
    when (not exists) $ error (printf "The file \"%s\" doesn't exists, when trying to read it" name')
    name <- newCString name'
    dataset <- newCString dataset'
    --Get the rank and dimensions, and calculate how many elements the dataset has
    rank <- getRankDataset name' dataset'
    dims <- getDimsDataset name' dataset' rank
    let total = SV.product dims
    -- Allocate enough memory for the array
    dat <- mallocArray (fromIntegral total) :: IO (Ptr b)
    -- Get the data via the FFI binding
    reader name dataset dat
    -- Make it to a foreign pointer, this will make sure it gets freed eventually and we can get a Vector out of it
    foreign_dat <- newForeignPtr finalizerFree (castPtr dat) :: IO (A.ForeignPtrs e)
    -- Give us back a vector :)
    rank <- getRankDataset name' dataset'
    dims <- getDimsDataset name' dataset' rank
    let askedRank = A.arrayRank (undefined :: A.Array sh b)
        shape = assert (if (askedRank == fromIntegral rank) then True else error (printf "Assertion failed for dataset %s in file %s. Expected rank: %i, actual rank: %i" dataset' name' askedRank (fromIntegral rank :: Int)))
                . A.listToShape . reverse . SV.toList . SV.map fromIntegral $ dims
    return $ A.fromForeignPtrs shape foreign_dat

readDatasets :: forall e a b sh.(Storable b, A.EltRepr b ~ e, A.Elt b, A.ForeignPtrs e ~ ForeignPtr a, A.Shape sh) =>
             (CString -> Ptr CString -> Ptr b -> IO ())
            -> String -> [String] -> IO (A.Array sh b)
readDatasets reader name' datasets' = do
    exists <- doesFileExist name'
    when (not exists) $ error (printf "The file \"%s\" doesn't exists, when trying to read it" name') 
    name <- newCString name'
    datasets <- mapM newCString datasets'
    datasetsPtr <- newArray0 nullPtr datasets
    --Get the rank and dimensions, and calculate how many elements the dataset has
    rank <- getRankDataset name' (head datasets')
    dims <- getDimsDataset name' (head datasets') rank
    let askedRank = A.arrayRank (undefined :: A.Array sh b)
        n = length datasets'
        total = assert (askedRank == 1 + fromIntegral rank) $ n * fromIntegral (SV.product dims)
        newdims = reverse . (n:) . SV.toList . SV.map fromIntegral $ dims
        shape = A.listToShape newdims
    -- Allocate enough memory for the array
    dat <- mallocArray (fromIntegral total) :: IO (Ptr b)
    -- Get the data via the FFI binding
    reader name datasetsPtr dat
    -- Make it to a foreign pointer, this will make sure it gets freed eventually and we can get a Vector out of it
    foreign_dat <- newForeignPtr finalizerFree (castPtr dat) :: IO (A.ForeignPtrs e)
    -- Give us back a vector :)
    return $ A.fromForeignPtrs shape foreign_dat

unsafeCastDataSet :: (A.Shape sh, A.Elt e1,  A.Elt e2, A.ForeignPtrs (A.EltRepr e1) ~ ForeignPtr a0, A.ForeignPtrs (A.EltRepr e2) ~ ForeignPtr b0) => 
    A.Array sh e1 -> A.Array sh e2
unsafeCastDataSet a = A.fromForeignPtrs (A.arrayShape a) $ castForeignPtr (A.toForeignPtrs a)

readDatasetInt64 :: (A.Shape sh) => String -> String -> IO (A.Array sh Int64)
readDatasetInt64 n1 n2 = do
    m <- readDataset c_readDatasetInt64 n1 n2
    return $ unsafeCastDataSet m

readDatasetDouble :: (A.Shape sh) => String -> String -> IO (A.Array sh Double)
readDatasetDouble n1 n2 = do
    m <- readDataset c_readDatasetDouble n1 n2
    return $ unsafeCastDataSet m

readDatasetComplex :: (A.Shape sh) => String -> String -> IO (A.Array sh ComplexDouble)
readDatasetComplex n1 n2 = do
    m <- readDataset c_readDatasetComplex n1 n2
    return $ unsafeCastDataSet m

readDatasetsDouble :: (A.Shape sh) => String -> [String] -> IO (A.Array sh Double)
readDatasetsDouble n ds = do
    m <- readDatasets c_readDatasetsDouble n ds
    return $ unsafeCastDataSet m

readDatasetsComplex :: (A.Shape sh) => String -> [String] -> IO (A.Array sh ComplexDouble)
readDatasetsComplex n ds = do
    m <- readDatasets c_readDatasetsComplex n ds
    return $ unsafeCastDataSet m

createDatasetDouble :: (A.Shape sh) => String -> String -> A.Array sh Double -> IO ()
createDatasetDouble name dataset dat= do
    let dims = SV.map fromIntegral . SV.fromList . reverse . A.shapeToList . A.arrayShape $ dat
    let rank = fromIntegral . SV.length $ dims
    let vec = SV.unsafeCast . A.toVectors $ dat
    createDatasetDoubleV name dataset rank dims vec

createDatasetComplex :: (A.Shape sh) => String -> String -> A.Array sh ComplexDouble -> IO ()
createDatasetComplex name dataset dat = do
    let dims = SV.map fromIntegral . SV.fromList . reverse . A.shapeToList . A.arrayShape $ dat :: SV.Vector CInt
    let rank = fromIntegral . SV.length $ dims
    let vec = SV.unsafeCast . A.toVectors $ dat
    createDatasetComplexV name dataset rank dims vec