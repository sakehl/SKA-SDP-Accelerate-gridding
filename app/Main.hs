module Main where

import ImageDataset
import System.Directory
import Control.Exception
import System.IO.Error hiding (catch)
import System.Environment
import Data.Array.Accelerate.LLVM.Native                  as CPU

data Args = Args { n   :: Maybe Int
                 , gpu :: Bool}

defArgs = Args (Just 500) False

main :: IO ()
main = do
    remover "result.h5"
    args <- getArgs
    let parsedArgs = parser defArgs args
        n_          = n parsedArgs
    fourier <- aw_gridding CPU.run CPU.runN "data/SKA1_Low_wkern2.h5" "data/SKA1_Low_akern3.h5" "data/SKA1_Low_quick.h5" n_
    putStrLn (show fourier)


remover :: FilePath -> IO ()
remover fileName = removeFile fileName `catch` handle
    where handle e
            | isDoesNotExistError e = return ()
            | otherwise = throwIO e

parser :: Args -> [String] -> Args
parser args []          = args
parser args ("-g":xs)   = parser args{gpu = True} xs
parser args ("-G":xs)   = parser args{gpu = True} xs
parser args ("-gpu":xs) = parser args{gpu = True} xs
parser args ("-GPU":xs) = parser args{gpu = True} xs
parser args ("-Gpu":xs) = parser args{gpu = True} xs
parser args ("-n":n:xs) = parser args{n = Just $ read n} xs
parser args ("-all":xs) = parser args{n = Nothing} xs
parser args xs          = error $ "Error while parsing" ++ show xs