import os

# ! Not yet ready to do something else than .jpg
# Function to rename multiple files 
def rename_folder_files(folder_path): 

    for count, filename in enumerate(os.listdir(folder_path)): 
        dst ="Hostel" + str(count) + ".jpg"
        src = folder_path + filename 
        dst = folder_path + dst 

        # rename() function will 
        # rename all the files 
        os.rename(src, dst) 
