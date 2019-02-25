{-# LANGUAGE ForeignFunctionInterface #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}

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

type Complex a = (a, a)
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
    c_readDatasetComplex :: CString -> CString -> Ptr ComplexDouble -> IO ()

foreign import ccall unsafe "createDatasetDouble"
    c_createDatasetDouble :: CString -> CString -> CInt -> Ptr CInt -> Ptr CDouble -> IO ()

foreign import ccall unsafe "createDatasetComplex"
    c_createDatasetComplex :: CString -> CString -> CInt -> Ptr CInt -> Ptr ComplexDouble -> IO ()



createh5File :: String -> IO ()
createh5File name' = do 
    name <- newCString name'
    c_createh5File name

createDatasetDouble :: String -> String -> CInt -> SV.Vector CInt -> SV.Vector CDouble -> IO ()
createDatasetDouble = createDataset c_createDatasetDouble
createDatasetComplex :: String -> String -> CInt -> SV.Vector CInt -> SV.Vector ComplexDouble -> IO ()
createDatasetComplex = createDataset c_createDatasetComplex

readDatasetDouble :: String -> String -> IO (SV.Vector CDouble)
readDatasetDouble = readDataset c_readDatasetDouble
readDatasetComplex :: String -> String -> IO (SV.Vector ComplexDouble)
readDatasetComplex = readDataset c_readDatasetComplex

createDataset :: Storable a => (CString -> CString -> CInt -> Ptr CInt -> Ptr a -> IO ()) -> 
    String -> String -> CInt -> SV.Vector CInt -> SV.Vector a -> IO ()
createDataset c_createDataset name' dataset' rank dims dat = do
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

readDataset :: Storable a => (CString -> CString -> Ptr a -> IO ()) 
            -> String -> String -> IO (SV.Vector a)
readDataset reader name' dataset' = do
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
readDataset2 :: (A.ForeignPtrs (A.EltRepr e) , Storable e, A.Elt e) => (CString -> CString -> Ptr e -> IO ()) 
            -> String -> String -> IO (A.Array A.DIM1 e)
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
    foreign_dat <- newForeignPtr finalizerFree dat :: IO (ForeignPtr e)
    -- Give us back a vector :)
    --return $ SV.unsafeFromForeignPtr0 foreign_dat (fromIntegral total)
    let shape = A.Z A.:.  (fromIntegral total) :: (A.:.) A.Z Int
    return $ A.fromForeignPtrs shape foreign_dat
-}
readDataset2 :: (CString -> CString -> Ptr CDouble -> IO ()) 
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
    foreign_dat <- newForeignPtr finalizerFree (castPtr dat) :: IO (ForeignPtr Double)
    -- Give us back a vector :)
    --return $ SV.unsafeFromForeignPtr0 foreign_dat (fromIntegral total)
    let shape = A.Z A.:.  (fromIntegral total) :: (A.:.) A.Z Int
    return $ A.fromForeignPtrs shape foreign_dat

readDatasetDoubleA ::  String -> String -> IO (A.Array A.DIM1 CDouble)
readDatasetDoubleA name dataset = do
    vector <- readDatasetDouble name dataset
    let vector2 = SV.unsafeCast vector
    let shape = A.Z A.:. SV.length vector2
    return $ A.fromVectors shape vector2

instance (Storable a => Storable (a, a)) where
    sizeOf _ = 2 * (sizeOf (undefined :: a))
    alignment _ = alignment (undefined :: a)
    peek ptr = do
        let dptr = castPtr ptr
        r <- peek dptr
        i <- peekElemOff dptr 1
        return (r,i)
    poke ptr (r,i) = do
        let dptr = castPtr ptr
        poke dptr r
        pokeElemOff dptr 1 i
        return ()


testDims :: SV.Vector CInt
testDims = SV.fromList [3,3]

testData :: SV.Vector CDouble
testData = SV.generate 9 fromIntegral

testDataC :: SV.Vector ComplexDouble
testDataC = SV.generate 9 (\x -> let x' = fromIntegral x in (x', x'))