# -*- coding: utf-8 -*-

"""
Created on Tue Dec  8 14:59:11 2020
@author: Jean
"""

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import pandas as pd
import os
from skimage import filters
import cv2
from streamlit_drawable_canvas import st_canvas
import threading
import pickle
import copy
import exifread


def main():
    st.set_page_config(page_title="Image viewer", page_icon=":microscope:",layout="wide")
    
    state=State()
    
    with st.sidebar:
        state.file_dir=st.text_input('File directory',"E:/optorhoa/201210_RPE1_optoRhoa_RBDiRFP/")
        try:
            state.filename=file_selector(state.file_dir)
        except:
            state.filename=None
            st.write("No such directory or no .nd file in this directory")
    
    if state.resultsfile==None:
        state.resultsfile='./pdresults.pkl'
        resultspd=pd.DataFrame()
        resultspd.to_pickle(state.resultsfile)
        with open('./results.pkl', 'wb') as output:
            results=[]
            pickle.dump(results, output, pickle.HIGHEST_PROTOCOL)
    
    if not state.filename==None:
        look_images(state)
    

    state.new_exp=False

class State:
    resultsfile=None
    temppos=0
    
class WL:
    def __init__(self,name,step=1):
        self.name=name
        self.step=step
        
class Exp:
    def __init__(self,expname,wl=[],nbpos=1,nbtime=1,comments=[]):
        self.name=expname
        self.nbpos=nbpos
        self.nbtime=nbtime
        self.wl=wl
        self.nbwl=len(wl)
        self.comments=comments
        if self.nbtime==1:
            self.timestep=0
        else:
            maxwl_ind=min(list(range(self.nbwl)), key=lambda ind:self.wl[ind].step)
            with open(self.get_image_name(maxwl_ind,timepoint=1), 'rb') as opened:
                tags = exifread.process_file(opened)
                time_str=tags['Image DateTime'].values
                h, m, s = time_str.split(' ')[1].split(':')
                time1=int(h) * 3600 + int(m) * 60 + float(s)
            with open(self.get_image_name(maxwl_ind,timepoint=int((nbtime-1)/self.wl[maxwl_ind].step+1)), 'rb') as opened:
                tags = exifread.process_file(opened)
                time_str=tags['Image DateTime'].values
                h, m, s = time_str.split(' ')[1].split(':')
                time2=int(h) * 3600 + int(m) * 60 + float(s)
            self.timestep=(time2-time1)/self.nbtime
    
    def get_image_name(self,wl_ind,pos=1,timepoint=1):
        if self.nbtime==1:
            tpstring=''
        else:
            tpstring='_t'+str(timepoint)
        if self.nbpos==1:
            posstring=''
        else:
            posstring='_s'+str(pos)
        return self.name+'_w'+str(wl_ind+1)+self.wl[wl_ind].name+posstring+tpstring+'.tif'
    
    def get_first_image(self,wl_ind,pos=1,timepoint=''):
        timepoint=1
        return Image.open(self.get_image_name(wl_ind,pos,timepoint))
    
    def get_last_image(self,wl_ind,pos=1,timepoint=1):
        last_ind=int(self.nbtime/self.wl[wl_ind].step-1)*self.wl[wl_ind].step+1
        return Image.open(self.get_image_name(wl_ind,pos,last_ind))

def get_exp(filename):
    nb_pos=1
    nb_wl=1
    with open(filename,'r') as file:
        i=0
        line=file.readline()
        comments=[]
        iscomments=False
        while not line.rstrip().split(', ')[0]=='"NTimePoints"' and i<50:
            if line.rstrip().split(', ')[0]=='"StartTime1"':
                iscomments=False
            if iscomments:
                comments.append(line.rstrip())
            if line.rstrip().split(', ')[0]=='"Description"':
                iscomments=True
                comments.append(str(line.rstrip().split(', ')[1]))
            line=file.readline()
            i+=1
        #get number of timepoints
        nb_tp=int(line.rstrip().split(', ')[1])
        line=file.readline()
        
        #get positions if exist
        if line.split(', ')[1].rstrip('\n')=='TRUE':
            line=file.readline()
            nb_pos=int(line.split(', ')[1].rstrip('\n'))
            for i in range(nb_pos):
                file.readline()            
        file.readline()
        
        #get number of wavelengths
        line=file.readline()
        nb_wl=int(line.rstrip().split(', ')[1])
    
        #create all new wavelengths
        wl=[]
        for i in range (nb_wl):
            line=file.readline()
            wl.append(WL(line.rstrip().split(', ')[1].strip('\"')))
            file.readline()
    
        #change the time steps
        line=file.readline()
        while line.split(', ')[0].strip('\"')=='WavePointsCollected':
            sep=line.rstrip().split(', ')
            if len(sep)>3:
                wl[int(sep[1])-1].step=int(sep[3])-int(sep[2])
            line=file.readline()
        
        expname=filename.rstrip('.nd')
        
        return Exp(expname,wl,nb_pos,nb_tp,comments)              

def look_images(state):
    with st.sidebar:
        exp=get_exp(state.filename)
        state.exp=exp
        st.markdown('**'+'Comments :'+'**')
        for comment in exp.comments:
            if not comment=='':
                st.markdown('_'+comment+'_') 
        st.write('**'+"Number of positions : "+'**'+str(exp.nbpos))
        st.write('**'+"Number of time steps : "+'**'+str(exp.nbtime))
        st.write('**'+"Time step : "+'**'+str(exp.timestep)+' sec')
        state.pos=st.selectbox('Position',range(1,exp.nbpos+1))
            
    if (not state.temppos==state.pos) or (not state.exp.name==state.tempexpname):
        state.new_exp=True
        state.temppos=state.pos
        state.tempexpname=state.exp.name
    
    #decide how many columns should be done       
    try:
        col=list(st.beta_columns(int((len(exp.wl)+1)/2)+1))
        state.wlthres=[None]*4
        state.fig=[None]*4 
        state.img=[None]*4
    except:
        st.write('Please load a valid experiment first')
        col=None
    
    #Create columns only if experiment is loaded
    if col is not None:
        
        for i in range(len(exp.wl)):
                
            #checking if I should update the temporary image
            
            if state.new_exp:
                create_image(state,i)
                #storing the image to be opened faster afterwards
                state.img[i]=Image.open('temp_wl_'+str(i)+'.png')
            else:
                if (not state.tempcoeff_seg==state.coeff_seg) and (i==state.wl_seg):
                    create_image(state,i)
                    state.tempcoeff_seg=state.coeff_seg
                    #storing the image to be opened faster afterwards
                    state.img[i]=Image.open('temp_wl_'+str(i)+'.png')
                if (not state.tempcoeff_act==state.coeff_act) and (i==state.wl_act) and (not state.draw):
                    create_image(state,i)
                    state.tempcoeff_act=state.coeff_act
                    #storing the image to be opened faster afterwards
                    state.img[i]=Image.open('temp_wl_'+str(i)+'.png')         
            
            
            #displaying intermediate columns
            with col[int(i/2)]:
                st.write(state.exp.wl[i].name)
                st.image(state.img[i], use_column_width=True)
    
        #last column
        with col[-1]:
               last_col(state)


def create_image(state,i):
    img1=np.array(state.exp.get_first_image(i,state.pos))
    img2=np.array(state.exp.get_last_image(i,state.pos))
        
    if 1==0:
        pass
    else:
        fig1 = plt.figure()
        plt.subplot(1,2,1)
        a=plt.imshow(img1,cmap='gray')
        a.axes.axis('off')
        plt.subplot(1,2,2)
        b=plt.imshow(img2,cmap='gray')
        b.axes.axis('off')
        plt.tight_layout()
    fig1.savefig('temp_wl_'+str(i)+'.png',bbox_inches="tight")

def last_col(state):
    inds=range(len(state.exp.wl))
    state.wls_mov=st.multiselect('Channels for movies',inds,format_func=lambda i: state.exp.wl[i].name,key='meas')
    if st.button('Make movies'):
        th=threading.Thread(target=make_all_movies,args=[state]) 
        th.start()

def make_all_movies(state):
    exp=copy.deepcopy(state.exp)
    if exp.nbpos==1:
        posstring=''
    else:
        posstring='_s'+str(state.pos)
    for i in state.wls_mov:
        make_movie(exp.name+'_w'+str(i+1)+exp.wl[i].name+posstring+'_t','.tif',exp.nbtime,1)

def make_movie(file1,file2,nbimg,step,crop=False):
    #initial step to find size
    img=Image.open(file1+str(1+2*step)+file2)
    size = img.size
    if crop:
        size=(1024,1024)
        start=512
        end=1537
        imgarray=np.array(img)[start:end,start:end]
    else:
        imgarray=np.array(img)
    sort=np.sort(imgarray.flatten())
    #sortcut=sort[int(len(sort)/1000):int(999*len(sort)/1000)]
    maxi,mini=np.max(sort),np.min(sort)
    
    out = cv2.VideoWriter(file1+'all.avi',0,7, size)
    for i in range(1,nbimg+1,step):
        if crop:
            img=np.array(Image.open(file1+str(i)+file2))[start:end,start:end]
        else:
            img=np.array(Image.open(file1+str(i)+file2))
        Image.fromarray(improve_contrast(img,mini,maxi)).save('./temp.png')
        img=cv2.imread('./temp.png')
        out.write(img)

def improve_contrast(img,mini,maxi):
    eq=255*(img.astype('float')-mini)/(maxi-mini)
    eq[eq<0]=0
    eq[eq>255]=255
    return np.uint8(eq)

def file_selector(folder_path='.',extension='.nd'):
    filenames = [f for f in os.listdir(folder_path) if f.endswith(extension)]
    selected_filename = st.selectbox('Select a file', filenames)
    return os.path.join(folder_path, selected_filename)            

    
if __name__ == "__main__":
   main()
   