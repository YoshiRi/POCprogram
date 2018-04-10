import cv2
import numpy as np
import imregpoc
import math

def ChangeContrast(image,a):
  lut = [ np.uint8(255.0 / (1 + math.exp(-a * (i - 128.) / 255.))) for i in range(256)] 
  result_image = np.array( [ lut[value] for value in image.flat], dtype=np.uint8 )
  result_image = result_image.reshape(image.shape)
  return result_image

def ChangeLight(image,a=1.2,b=15):
  lut = [ max(0,min(255,(np.uint8( a*i+b )))) for i in range(256)] 
  result_image = np.array( [ lut[value] for value in image.flat], dtype=np.uint8 )
  result_image = result_image.reshape(image.shape)
  return result_image


from scipy import ndimage

#A=ndimage.rotate(cv2.imread('clouds13.png',0),270,reshape=False)
A=np.rot90(np.rot90(np.rot90(cv2.imread('clouds13.png',0))))
B=cv2.imread('clouds37.png',0)

A1=ChangeContrast(A,15)
A2=ChangeLight(A)

poc=imregpoc.imregpoc(A,B)
poc.stitching()


poc1=imregpoc.imregpoc(A1,B)
poc1.stitching()

poc2=imregpoc.imregpoc(A2,B)
poc2.stitching()

fp=imregpoc.TempMatcher(A,'SIFT')
fp.match(B)
fp.stitching()


fp1=imregpoc.TempMatcher(A1,'SIFT')
fp1.match(B)
fp1.stitching()

fp2=imregpoc.TempMatcher(A2,'SIFT')
fp2.match(B)
fp2.stitching()






