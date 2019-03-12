module Main where
    
import ImageDataset
import System.Directory
import Control.Exception
import System.IO.Error hiding (catch)
import System.Environment

main :: IO ()
main = do 
    remover "result.h5"
    args <- getArgs
    let n = case args of
                [] -> Just 50
                "all":_ -> Nothing
                x:_ -> Just $ read x
    fourier <- aw_gridding "data/SKA1_Low_wkern2.h5" "data/SKA1_Low_akern2.h5" "data/SKA1_Low_quick.h5" n
    putStrLn (show fourier)


remover :: FilePath -> IO ()
remover fileName = removeFile fileName `catch` handle
    where handle e
            | isDoesNotExistError e = return ()
            | otherwise = throwIO e