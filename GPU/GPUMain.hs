module Main where
    
import ImageDataset
import System.Directory
import Control.Exception
import System.IO.Error hiding (catch)
import System.Environment
import Data.Array.Accelerate.LLVM.Native                  as CPU
import Data.Array.Accelerate.LLVM.PTX                     as GPU


import Data.Array.Accelerate.Debug                        as A

main :: IO ()
main = do 
    remover "result.h5"
    --setFlag dump_phases
    args <- getArgs
    let n = case args of
                [] -> Just 500
                "all":_ -> Nothing
                x:_ -> Just $ read x
    let isGPU = case args of
        [] -> False
        [_:"GPU":_] -> True
        [_:"Gpu":_] -> True
        [_:"gpu":_] -> True
        _           -> False
    fourier <- if isGPU 
                then aw_gridding GPU.run GPU.runN "data/SKA1_Low_wkern2.h5" "data/SKA1_Low_akern3.h5" "data/SKA1_Low_quick.h5" n
                else aw_gridding CPU.run CPU.runN "data/SKA1_Low_wkern2.h5" "data/SKA1_Low_akern3.h5" "data/SKA1_Low_quick.h5" n
    putStrLn (show fourier)


remover :: FilePath -> IO ()
remover fileName = removeFile fileName `catch` handle
    where handle e
            | isDoesNotExistError e = return ()
            | otherwise = throwIO e