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

import qualified Data.Vector.Storable as SV
import qualified Data.Array.Accelerate.IO.Data.Vector.Storable as A
import qualified Data.Array.Accelerate                         as A
import qualified Data.Array.Accelerate.IO.Foreign.ForeignPtr   as A
import qualified Data.Array.Accelerate.Array.Sugar             as A hiding(shape)
import Data.Array.Accelerate.Data.Complex

import Debug.Trace
--type Complex a = (a, a)
type ComplexCDouble = Complex CDouble
type ComplexDouble = Complex Double

---- Actual functionallity
foreign import ccall unsafe "createh5File"
    c_createh5File :: CString -> IO ()

foreign import ccall unsafe "getRankDataset"
    c_getRankDataset :: CString -> CString -> IO CInt

foreign import ccall unsafe "getDimsDataset"
    c_getDimsDataset :: CString -> CString -> CInt -> Ptr CInt -> IO ()

foreign import ccall unsafe "readDatasetDouble"
    c_readDatasetDouble :: CString -> CString -> Ptr CDouble -> IO ()

foreign import ccall unsafe "readDatasetComplex"
    c_readDatasetComplex :: CString -> CString -> Ptr ComplexCDouble -> IO ()

foreign import ccall unsafe "createDatasetDouble"
    c_createDatasetDouble :: CString -> CString -> CInt -> Ptr CInt -> Ptr CDouble -> IO ()

foreign import ccall unsafe "createDatasetComplex"
    c_createDatasetComplex :: CString -> CString -> CInt -> Ptr CInt -> Ptr ComplexCDouble -> IO ()



createh5File :: String -> IO ()
createh5File name' = do 
    name <- newCString name'
    c_createh5File name

createDatasetDoubleV :: String -> String -> CInt -> SV.Vector CInt -> SV.Vector CDouble -> IO ()
createDatasetDoubleV = createDatasetV c_createDatasetDouble
createDatasetComplexV :: String -> String -> CInt -> SV.Vector CInt -> SV.Vector ComplexCDouble -> IO ()
createDatasetComplexV = createDatasetV c_createDatasetComplex

readDatasetDoubleV :: String -> String -> IO (SV.Vector CDouble)
readDatasetDoubleV = readDatasetV c_readDatasetDouble
readDatasetComplexV :: String -> String -> IO (SV.Vector ComplexCDouble)
readDatasetComplexV = readDatasetV c_readDatasetComplex

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

readDatasetV :: Storable a => (CString -> CString -> Ptr a -> IO ()) 
            -> String -> String -> IO (SV.Vector a)
readDatasetV reader name' dataset' = do
    name <- newCString name'
    dataset <- newCString dataset'
    --Get the rank and dimensions, and calculate how many elements the dataset has
    rank <- getRankDataset name' dataset'
    dims <- getDimsDataset name' dataset' rank
    let total = SV.product dims
    -- Allocate enough memory for the array
    dat <- mallocArray (fromIntegral total)
    -- Get the data via the FFI binding
    reader name dataset dat
    -- Make it to a foreign pointer, this will make sure it gets freed eventually and we can get a Vector out of it
    foreign_dat <- newForeignPtr finalizerFree dat
    -- Give us back a vector :)
    return $ SV.unsafeFromForeignPtr0 foreign_dat (fromIntegral total)

{-
readDataset3 :: (Storable e, A.Elt e) =>
             (CString -> CString -> Ptr e -> IO ())
            -> String -> String -> IO (A.Array A.DIM1 e)
readDataset3 reader name' dataset' = do
    name <- newCString name'
    dataset <- newCString dataset'
    --Get the rank and dimensions, and calculate how many elements the dataset has
    rank <- getRankDataset name' dataset'
    dims <- getDimsDataset name' dataset' rank
    let total = SV.product dims
    -- Allocate enough memory for the array
    dat <- mallocArray (fromIntegral total) :: IO (Ptr e)
    -- Get the data via the FFI binding
    reader name dataset dat
    -- Make it to a foreign pointer, this will make sure it gets freed eventually and we can get a Vector out of it
    foreign_dat <- newForeignPtr finalizerFree (castPtr dat) :: IO (A.ForeignPtrs (A.EltRepr e))
    -- Give us back a vector :)
    --return $ SV.unsafeFromForeignPtr0 foreign_dat (fromIntegral total)
    let shape = A.Z A.:.  (fromIntegral total) :: (A.:.) A.Z Int
    return $ A.fromForeignPtrs shape foreign_dat
-}

readDataset2 :: --(Storable e, A.Elt e, A.ForeignPtrs e ~ a) =>
             (CString -> CString -> Ptr CDouble -> IO ())
            -> String -> String -> IO (A.Array A.DIM1 CDouble)
readDataset2 reader name' dataset' = do
    name <- newCString name'
    dataset <- newCString dataset'
    --Get the rank and dimensions, and calculate how many elements the dataset has
    rank <- getRankDataset name' dataset'
    dims <- getDimsDataset name' dataset' rank
    let total = SV.product dims
    -- Allocate enough memory for the array
    dat <- mallocArray (fromIntegral total)
    -- Get the data via the FFI binding
    reader name dataset dat
    -- Make it to a foreign pointer, this will make sure it gets freed eventually and we can get a Vector out of it
    foreign_dat <- newForeignPtr finalizerFree (castPtr dat) :: IO (A.ForeignPtrs (A.EltRepr CDouble))
    -- Give us back a vector :)
    --return $ SV.unsafeFromForeignPtr0 foreign_dat (fromIntegral total)
    let shape = A.Z A.:.  (fromIntegral total)
    return $ A.fromForeignPtrs shape foreign_dat

readDatasetDouble ::  (A.Shape sh) => String -> String -> IO (A.Array sh Double)
readDatasetDouble name dataset = do
    vector <- readDatasetDoubleV name dataset
    let vectorD = SV.unsafeCast vector

    rank <- getRankDataset name dataset
    dims <- getDimsDataset name dataset rank
    let shape = A.listToShape . SV.toList . SV.map fromIntegral $ dims
    return $ A.fromVectors shape vectorD

readDatasetComplex :: (A.Shape sh) => String -> String -> IO (A.Array sh ComplexDouble)
readDatasetComplex name dataset = do
    vector <- readDatasetComplexV name dataset
    let vectorD = SV.unsafeCast vector

    rank <- getRankDataset name dataset
    dims <- getDimsDataset name dataset rank
    let shape = A.listToShape . SV.toList . SV.map fromIntegral $ dims
    return $ A.fromVectors shape vectorD

createDatasetDouble :: (A.Shape sh) => String -> String -> A.Array sh Double -> IO ()
createDatasetDouble name dataset dat= do
    let dims = SV.map fromIntegral . SV.fromList . A.shapeToList . A.arrayShape $ dat
    let rank = fromIntegral . SV.length $ dims
    let vec = SV.unsafeCast . A.toVectors $ dat
    createDatasetDoubleV name dataset rank dims vec

createDatasetComplex :: (A.Shape sh) => String -> String -> A.Array sh ComplexDouble -> IO ()
createDatasetComplex name dataset dat = do
    let dims = SV.map fromIntegral . SV.fromList . A.shapeToList . A.arrayShape $ dat :: SV.Vector CInt
    let rank = fromIntegral . SV.length $ dims
    let vec = SV.unsafeCast . A.toVectors $ dat
    createDatasetComplexV name dataset rank dims vec