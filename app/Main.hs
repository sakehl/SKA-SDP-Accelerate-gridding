module Main where
    
import ImageDataset

main :: IO ()
main = do 
    fourier <- aw_gridding "data/SKA1_Low_wkern.h5" "data/SKA1_Low_akern.h5" "data/SKA1_Low_quick.h5"
    putStrLn (show fourier)