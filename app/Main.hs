{-# LANGUAGE CPP #-}
module Main where

import ImageDataset
import System.Directory
import Control.Exception
import System.IO.Error hiding (catch)
import System.Environment
import qualified Data.Array.Accelerate.LLVM.Native      as CPU
import qualified Data.Array.Accelerate.Interpreter      as I

#ifdef ACCELERATE_LLVM_PTX_BACKEND
import qualified Data.Array.Accelerate.LLVM.PTX         as GPU
#endif

data Args = Args { n   :: Maybe Int
                 , runner :: Run
                 , input :: String
                 , out :: Maybe String}

data Run = GPU | CPU | Inter

defArgs = Args (Just 1) CPU "data" Nothing

main :: IO ()
main = do
    --remover "result.h5"
    args <- getArgs
    let parsedArgs = parser defArgs args
        n_         = n parsedArgs
        backend    = runner parsedArgs
        dir        = input parsedArgs
        wkern      = dir ++ "/" ++ "SKA1_Low_wkern2.h5"
        akern      = dir ++ "/" ++ "SKA1_Low_akern3.h5"
        visdata    = dir ++ "/" ++ "SKA1_Low_quick.h5"
        outFile    = out parsedArgs
    case outFile of
        Nothing -> return ()
        Just fn  -> remover fn
    fourier <- case backend of
        CPU -> aw_gridding CPU.run CPU.runN wkern akern visdata n_ outFile
        Inter -> aw_gridding I.run I.runN wkern akern visdata n_ outFile
#ifdef ACCELERATE_LLVM_PTX_BACKEND
        GPU -> aw_gridding GPU.run GPU.runN wkern akern visdata n_ outFile
#endif
        _ -> aw_gridding CPU.run CPU.runN wkern akern visdata n_ outFile
    putStrLn (show fourier)


remover :: FilePath -> IO ()
remover fileName = removeFile fileName `catch` handle
    where handle e
            | isDoesNotExistError e = return ()
            | otherwise = throwIO e

parser :: Args -> [String] -> Args
parser args []            = args
parser args ("-debug":xs) = parser args{runner = Inter} xs
parser args ("-g":xs)     = parser args{runner = GPU} xs
parser args ("-G":xs)     = parser args{runner = GPU} xs
parser args ("-gpu":xs)   = parser args{runner = GPU} xs
parser args ("-GPU":xs)   = parser args{runner = GPU} xs
parser args ("-Gpu":xs)   = parser args{runner = GPU} xs
parser args ("-n":n:xs)   = parser args{n = Just $ read n} xs
parser args ("-all":xs)   = parser args{n = Nothing} xs
parser args ("-i":inp:xs) = parser args{input = inp} xs
parser args ("-o":out:xs) = parser args{out = Just out} xs
parser args xs            = error $ "Error while parsing" ++ show xs