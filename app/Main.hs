module Main where
    
import ImageDataset
import System.Directory
import Control.Exception
import System.IO.Error hiding (catch)
import System.Environment
import Data.Array.Accelerate.LLVM.Native                  as CPU

main :: IO ()
main = do 
    remover "result.h5"
    args <- getArgs
    let n = case args of
                [] -> Just 500
                "all":_ -> Nothing
                x:_ -> Just $ read x
    fourier <- aw_gridding CPU.run CPU.runN "data/SKA1_Low_wkern2.h5" "data/SKA1_Low_akern3.h5" "data/SKA1_Low_quick.h5" n 
    putStrLn (show fourier)


remover :: FilePath -> IO ()
remover fileName = removeFile fileName `catch` handle
    where handle e
            | isDoesNotExistError e = return ()
            | otherwise = throwIO e