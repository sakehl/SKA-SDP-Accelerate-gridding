module Main where
    
import ImageDataset
import System.Directory
import Control.Exception
import System.IO.Error hiding (catch)
import System.Environment
import Data.Array.Accelerate.LLVM.Native                  as CPU
import Data.Array.Accelerate.LLVM.PTX                     as GPU


import Data.Array.Accelerate.Debug                        as A

data Args = Args { n   :: Maybe Int
                 , gpu :: Bool
                 , input :: String}

defArgs = Args (Just 500) False "data"

main :: IO ()
main = do 
    remover "result.h5"
    --setFlag dump_phases
    args <- getArgs
    let parsedArgs = parser defArgs args
        n_         = n parsedArgs
        isGPU      = gpu parsedArgs
        dir        = input parsedArgs
        wkern      = dir ++ "/" ++ "SKA1_Low_wkern2.h5"
        akern      = dir ++ "/" ++ "SKA1_Low_akern3.h5"
        visdata    = dir ++ "/" ++ "SKA1_Low_quick.h5"
    fourier <- if isGPU 
                then aw_gridding GPU.run GPU.runN wkern akern visdata n_
                else aw_gridding CPU.run CPU.runN wkern akern visdata n_
    putStrLn (show fourier)


remover :: FilePath -> IO ()
remover fileName = removeFile fileName `catch` handle
    where handle e
            | isDoesNotExistError e = return ()
            | otherwise = throwIO e

parser :: Args -> [String] -> Args
parser args []            = args
parser args ("-g":xs)     = parser args{gpu = True} xs
parser args ("-G":xs)     = parser args{gpu = True} xs
parser args ("-gpu":xs)   = parser args{gpu = True} xs
parser args ("-GPU":xs)   = parser args{gpu = True} xs
parser args ("-Gpu":xs)   = parser args{gpu = True} xs
parser args ("-n":n:xs)   = parser args{n = Just $ read n} xs
parser args ("-all":xs)   = parser args{n = Nothing} xs
parser args ("-i":inp:xs) = parser args{input = inp} xs
parser args xs            = error $ "Error while parsing" ++ show xs