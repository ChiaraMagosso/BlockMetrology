## 
## Copyright (C) 2024 by
## Chiara Magosso
## 
## This work is licensed under a  
## 
## Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
## ( http://creativecommons.org/licenses/by-nc-sa/4.0/ )
## 
## Please contact chiara.magosso@polito.it for information.
##

import tifffile
import os
import codecs

def metadata(sem_image, path, name):
    metadata = open(f'{path}/{name}_metadata.txt', 'w')
    with tifffile.TiffFile(sem_image) as tif:
        for page in tif.pages:
            for tag in page.tags:
                tag_name, tag_value = tag.name, tag.value
                if  isinstance(tag_value, dict):
                    metadata.write(f'{tag_name} \n') 
                    for key, value in tag_value.items():
                        metadata.write(f'{key} \t {value} \n')
                else:
                    metadata.write(f"{tag_name} \t {codecs.decode(tag_value, 'unicode_escape') if isinstance(tag_value, str) else tag_value} \n")
    metadata.close()


path = 'C:/Users/machinelearning/Desktop/gui_metadata/' #modify accordingly 
for file in os.listdir(path):
    if file.endswith(".tif"):

        name = os.path.splitext(file)[0]
        print('\n'f'{name}')
        sem_image = (f'{path}{file}')
        
        print('Metadata extraction') 
        metadata(sem_image, path, name) 