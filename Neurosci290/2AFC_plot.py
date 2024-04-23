# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 17:12:28 2016

@author: dylantamayo
"""

import Tkinter
import tkFileDialog
import os, os.path
import matplotlib.pyplot as plt


def chooseFile():
    Tkinter.Tk().withdraw() # Close the root window
    thePath = tkFileDialog.askopenfilename()
    return thePath

def chooseDirectory():
    Tkinter.Tk().withdraw() # Close the root window
    thePath = tkFileDialog.askdirectory(title='Please select a directory')
  
    return thePath

def processFiles():
    dir = chooseDirectory()
#    dir = '/Users/dylantamayo/Desktop/Independent Project Stuff'
    allFiles = []
    for dirName, subList, fileList in os.walk(dir):
        for f in fileList:
            if f.startswith("."):
                continue
            fullname = os.path.join(dirName,f)
            if f.endswith('.csv'):
                allFiles.append(fullname)
    print "Number of Subjects:",len(allFiles)
    return allFiles

def listcompiler(files):
    biglist = []
    list1_10 = []
    list11_20 = []
    list21_30 = []
    list31_40 = []
    list41_50 = []
    for onefile in files:
        csvfile = open(onefile, 'rb')
        rownum = 1
        for row in csvfile:
            datalist = row.split(",")
            try:
                val1 = float((datalist[7]))
                biglist.append((val1,int(datalist[2])))
                if rownum < 11:
                    list1_10.append(int(datalist[2]))
                if rownum >= 11 and rownum < 21:
                    list11_20.append(int(datalist[2]))
                if rownum >= 21 and rownum < 31:
                    list21_30.append(int(datalist[2]))
                if rownum >= 31 and rownum < 41:
                    list31_40.append(int(datalist[2]))
                if rownum >= 41 and rownum < 51:
                    list41_50.append(int(datalist[2]))
                rownum += 1
                
            except ValueError:
                continue
    return biglist,list1_10,list11_20,list21_30,list31_40,list41_50

def dataplotter(files):
    biglist,a,b,c,d,e = listcompiler(files)
    xval = [1.25,2.5,5.0,10.0,20.0]
    yval1 = [y[1] for y in biglist if y[0] == 1.25]
    yval2 = [y[1] for y in biglist if y[0] == 2.5]
    yval3 = [y[1] for y in biglist if y[0] == 5.0]
    yval4 = [y[1] for y in biglist if y[0] == 10.0]
    yval5 = [y[1] for y in biglist if y[0] == 20.0]
    percentdiff = [sum(yval1)*100.0/len(yval1),sum(yval2)*100.0/len(yval2),
                   sum(yval3)*100.0/len(yval3),sum(yval4)*100.0/len(yval4),
                   sum(yval5)*100.0/len(yval5)]
    
    totalcorrect = (sum(yval1)+sum(yval2)+sum(yval3)+
                    sum(yval4)+sum(yval5))*100.0/len(biglist)
    print "Total percent correct:",totalcorrect
    print percentdiff
    plt.title('''Accuracy in Distinguishing Between
    Two Shades of Blue''')
    plt.xlabel('Percent Difference (%)')
    plt.ylabel('Percent Correct (%)')
    plt.plot(xval,percentdiff,'bo-')
    plt.gca().invert_xaxis()
    plt.show()
    plotter2(a,b,c,d,e)

def plotter2(a,b,c,d,e):
    yval1 =[sum(a)*100.0/len(a),sum(b)*100.0/len(b),sum(c)*100.0/len(c),
            sum(d)*100.0/len(d),sum(e)*100.0/len(e)]
    xval = [10,20,30,40,50]
    plt.title('''Total Accuracy as Experiment Progresses''')
    plt.xlabel('Trial Range')
    plt.ylabel('Percent Correct (%)')
    plt.plot(xval,yval1,'ro-')
    plt.show()


if __name__ == "__main__":
    files = processFiles()
    dataplotter(files)