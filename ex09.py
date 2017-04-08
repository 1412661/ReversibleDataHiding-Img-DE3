#!/home/quocbao/anaconda3/bin/python3.5
import numpy as np
import bitarray
from scipy.misc import imread, imsave
import matplotlib.pyplot as plt
import zlib

def isExpandable(l,h):
    if abs(2*h)   <= min(2*(255-l), 2*l+1) and \
       abs(2*h+1) <= min(2*(255-l), 2*l+1):
            return 1
    return 0

def isChangeable(l,h):
    if abs(2*(h//2))   <= min(2*(255-l), 2*l+1) and \
       abs(2*(h//2)+1) <= min(2*(255-l), 2*l+1):
            return 1
    return 0
    

def embed_DE(cover_img_file, msg_file, stego_img_file):
    """
    Embeds a message into a cover image using the DE (Difference Expansion) method.
    
    Parameters
    ----------
    cover_img_file : string
        The name of the cover img file.
    msg_file : string
        The name of the message file.
    stego_img_file : string
        The name of the stego img file.
    
    Returns
    -------
    result : bool
        result = True if it's possible to embed; False otherwise.
    """
    # TODO
        
    cover_img = imread(cover_img_file)
    height = cover_img.shape[0]
    width  = cover_img.shape[1]
    
    # Seperate original image into pairs of pixel
    cover_img = cover_img.reshape(cover_img.size)
    
    f = open(msg_file, 'r')
    msg = f.read()
    f.close()
    
    msg_bits = bitarray.bitarray()
    msg_bits.frombytes(msg.encode())
    
    de = list()
    ez = list()
    en = list()
    cn = list()
    nc = list()
    
    for i in range(0, cover_img.size, 2):
        l = (int(cover_img[i]) + int(cover_img[i+1])) // 2
        h = int(cover_img[i]) - int(cover_img[i+1])

        # Step 1: Calculate l and h
        de.append([l,h])
        
        # Step 2: Create four disjoint sets of difference values
        if isExpandable(l,h):
            if h in range(-1,1):        # h = 0 or h = -1
                ez.append(i//2)         # it mean [i] and [i+1]
            else:
                en.append(i//2)
        else:
            if isChangeable(l,h):
                cn.append(i//2)
            else:
                nc.append(i//2)
                
    print("Size EZ: ",len(ez))
    print("Size EN: ",len(en))
    print("Size CN: ",len(cn))
    print("Size NC: ",len(nc))
    
    # Check capacity
    print("Size of message: ", len(msg_bits)," bit(s)")
    if len(msg_bits) > len(ez+en+cn):
        print("Can't embed the message !")
        return 0
    else:
        print("We can embed the message !")
        
    # Step 3: Create location map
    lm = bitarray.bitarray('0' * (cover_img.size//2))

    # Mark location map
    for p in ez+en:
        lm[p] = '1'
    t = lm
    
    # Compress location map      
    zip_lm = bitarray.bitarray()
    zip_lm.frombytes(zlib.compress(lm.tobytes()))
    print("Location map before compression: ", len(lm))
    print("Location map after compression:  ", len(zip_lm))


    # Step 4: Collect original LSBs
    lsb = bitarray.bitarray()
    for i in cn:
        h = de[i][1] 
        if  h != -2 and h != 1:
            lsb.extend(str(de[i][1] & 1))
    print("Number of h = 1 or h = -2 in CN: ",len(cn) - len(lsb))
    
    # Step 5: Add padding to payload
    msg_bits.extend('1' + '0'*(len(ez+en+cn) - len(zip_lm) - len(lsb) - len(msg_bits) - 1))
        
    #print(len(msg_bits) + len(zip_lm) + len(lsb))
    #print(len(ez+en+cn))
    
    # Step 6: Create bitstream
    bitstream = zip_lm + lsb + msg_bits
    c = 0
    
    # Step 7: Embed bitstream into cover image
    stego_img = cover_img
    for i in range(0, stego_img.size, 2):        
        l = (int(stego_img[i]) + int(stego_img[i+1]))//2
        h = int(stego_img[i]) - int(stego_img[i+1])
                
        if isExpandable(l,h):
            h = int(2*h) + int(bitstream[c])
            c = c + 1
        elif isChangeable(l,h):
            h = int(2*(h//2)) + int(bitstream[c])
            c = c + 1
                
        stego_img[i]   = l + int((h+1)//2)
        stego_img[i+1] = l - int(h//2)

    # Step 8: Save stego image
    stego_img = stego_img.reshape(height,width)
    imsave(stego_img_file, stego_img)
    
    print("[Notice] Embeding process completed !\n")
        
        
def extract_DE(stego_img_file, extracted_msg_file, recovered_cover_img):
    """
    Extracts the message from a stego image using the DE method.
    
    Parameters
    ----------
    stego_img_file : string
        The name of the stego img file.
    extracted_msg_file : string
        The name of the extracted message file.
    recovered_cover_img_file : string
        The name of the recovered cover img file.
    """
    # TODO
    
    stego_img = imread(stego_img_file)
    height = stego_img.shape[0]
    width  = stego_img.shape[1]
    
    # Seperate original image into pairs of pixel
    stego_img = stego_img.reshape(stego_img.size)
    
    de = list()
    ch = list()
    nc = list()
    
    for i in range(0, stego_img.size, 2):
        l = (int(stego_img[i]) + int(stego_img[i+1])) // 2
        h = int(stego_img[i]) - int(stego_img[i+1])

        # Step 1: Calculate l and h and store it into de[]
        de.append([l,h])
        
        # Step 2: Create two disjoint sets of difference values
        if isChangeable(l,h):
            ch.append(i//2)
        else:
            nc.append(i//2)
        
    # Step 3: Collect bitstream
    bitstream = bitarray.bitarray()
    for i in ch:
        bitstream.extend(str(de[i][1] & 1))

    
    # Step 4: Decode location map
    unzip_lm = bitarray.bitarray()
    unzip_lm.frombytes(zlib.decompress(bitstream.tobytes()))

    # Compress location map again to get the size of compressed location map
    zip_lm = bitarray.bitarray()
    zip_lm.frombytes(zlib.compress(unzip_lm.tobytes()))
    
    # Get original LSBs and payload
    data = bitstream[len(zip_lm):]
    lsb_index = 0
    
    # Step 5: Restore original image
    recovered_img = stego_img
    for i in ch:
        l = de[i][0]
        h = de[i][1]
        
        if unzip_lm[i]:
            h = h // 2                  # DE
        else:
            if 0 <= h and h <= 1:
                h = 1
            elif -2 <= h and h <= -1:
                h = -2
            else:
                h = 2*(h//2) + data[lsb_index]  # LSB
                lsb_index += 1
        
        # Restore original image
        recovered_img[i*2]   = l + int((h+1)//2)
        recovered_img[i*2+1] = l - int(h//2)
    
    # Write down recovered image  
    recovered_img = recovered_img.reshape(height,width)
    imsave(recovered_cover_img, recovered_img)
         
    # Step 6: Get payload 
    payload = data[lsb_index:]                  # Cut original LSBs
    while not (payload[len(payload) - 1]):
        payload.pop()
    payload.pop()
    
    msg = payload.tostring()
    
    # Write message
    f = open(extracted_msg_file, 'w')
    f.write(msg)
    f.close()   
    
    print("[Notice] Extracting process completed !\n")

    

embed_DE("cover.bmp", "msg.txt", "stego.bmp")
extract_DE("stego.bmp", "hidden.txt", "recovered.bmp")
a = imread("cover.bmp")
b = imread("recovered.bmp")
if np.array_equal(a,b):
    print("[Notice] Recovered image look exactly as cover image")
else:
    print("[Notice] Recovered image is different from cover image")
    
f = open("msg.txt", 'r')
embeded_msg = f.read()
f.close()

f = open("hidden.txt", 'r')
extracted_msg = f.read()
f.close()

print("Embeded message:", embeded_msg)
print("Extracted message:", extracted_msg)
