# -*- coding: utf-8 -*-
"""
Created on Sat Apr 16 16:20:21 2016

@author: dylantamayo
"""
'''
Cleaned up some excess lines.
'''
from psychopy import visual, core, event, gui
import math, random, itertools

# GLOBAL state/variables
dataFile = None
mywin = None


# GLOBAL constants

TIME_THRESHOLD = 2    #trial stays "open" this long
BAD_NUMBER = 0.0
NUMBER_LOOPS = 5 #5 = 10 per percent difference

def getSettings():
    myDlg = gui.Dlg(title='Experiment')
    myDlg.addText('Subject info')
    myDlg.addField('Initials:')
    myDlg.addField('Gender:', choices=['M','F','O'])
#    myDlg.addField('Gender:', "M/F/O")
#    myDlg.addField('Experiment:', choices=["1", "2"])
    myDlg.show()  # show dialog and wait for OK or Cancel
    if myDlg.OK:  # then the user pressed OK
        thisInfo = myDlg.data
        return(tuple(thisInfo))
    else:
        print('user cancelled')

def initialize():
    global mywin#,myMouse
    global dataFile
    
    fileName = getSettings()##
    formatFile = '%s_%s' % (fileName)##   #'%s%s%s' % (fileName)
    dataFile = open(formatFile+'.csv', 'w')## #a simple text file with 'comma-separated-values'
    dataFile.write('Key pressed,Correct?,y-value,Red,Green,Left Blue,Right Blue,Percent Different,time\n')##
    
    mywin = visual.Window([800,600],monitor="testMonitor", units="deg", allowGUI=True, fullscr=True)

def inter_trial():
    x = random.randrange(2,6) #interval [n,x)
    betweentrials = visual.TextStim(mywin, color=1,wrapWidth=100,text="wait for next trial")
    betweentrials.draw()
    mywin.flip()
    core.wait(x/2) #CHANGE WAITING TIME BETWEEN TRIALS

def opening():
    display_message = """
    If the left square is lighter, press F.
    If the right square is lighter, press J.
              There are 50 trials.
       You have 2 seconds for each trial.
             Press Space to begin."""
    message = visual.TextStim(mywin,color=1, pos=(-0,0), wrapWidth=100,
                              text=display_message)
    message.draw()
    mywin.flip()
    event.waitKeys()

def logfunc(values, index):
    print values[index]
    return values[index]

def experiment(red,green,Lblue,Rblue): #logs the data for each trial
    x,time = expDisplay(red,green,Lblue,Rblue)
    gvalues = [('correct',1),('incorrect',0),('timeout',0)]
    lvalues = [('incorrect',0),('correct',1),('timeout',0)]
    if x == 'f': xindex = 0
    elif x == 'j': xindex = 1
    else: xindex = 2
    if Lblue > Rblue: #lighter
        y = logfunc(gvalues,xindex)
    else: #darker
        y = logfunc(lvalues,xindex)
    percentdiff = (abs(Lblue-Rblue))/.5*100
    print "logging",x,y[0],y[1],red,green,Lblue,Rblue,percentdiff,time
    datainfo = "%s,%s,%d,%1.3f,%1.3f,%1.3f,%1.3f,%1.3f,%1.3f\n" % (x,y[0],y[1],red,green,Lblue,Rblue,percentdiff,time) ##
    dataFile.write(datainfo) ##

def experiment1list(length):
    '''creates the list of trials, randomizes whether
    constant blue is left or right square'''
    x2 = [-.1,-.05,-.025,-.0125,-.00625,.00625,.0125,.025,.05,.1] #easier:[-.15,-.1,-.05,.05,.1,.15]
    y1 = []
    y2 = []
    for num in range(length):
        y1.append(x2)
        chain = list(itertools.chain.from_iterable(y1))
    for num in range(len(chain)):
        y2.append((chain[num],random.randrange(0,2)))
    print y2
    return y2

def experiment1(numtrials): #randomizes the trials, runs the experiment
    x = experiment1list(numtrials)
    random.shuffle(x)
    trialnum = 1
    for num in x:
        print num, trialnum
        if num[1] == 0:
            experiment(-.65,-.65,.5,.5+num[0])
        if num[1] == 1:
            experiment(-.65,-.65,.5+num[0],.5)
        trialnum += 1
        inter_trial()

def expDisplay(red,green,Lblue,Rblue):
    event.clearEvents()
    square1 = visual.Rect(win=mywin,width=7.5, height=7.5, fillColor=(red,green,Lblue), pos=[-10,0])
    square2 = visual.Rect(win=mywin,width=7.5, height=7.5, fillColor=(red,green,Rblue), pos=[10,0])
    '''towards -1 = darker, towards 1 = lighter. difference of .01 is somewhat distinguishable.
    '''
    clicked = False
    timeout = False
    start_time = core.getTime()
    
    while True: #this creates a never-ending loop

        square1.draw()
        square2.draw()
        flip_time = mywin.flip()
        
        if flip_time - start_time > TIME_THRESHOLD:
            timeout = True
            break

        keyboardpress = event.getKeys(keyList=['f','j'])

        if keyboardpress:
            clicked = True
            break
     
    time = flip_time - start_time
    if clicked:
        key = keyboardpress[0]
        return (key,time)
    else:
        return(BAD_NUMBER,time)

if __name__ == "__main__":
    initialize()
    opening()
    experiment1(NUMBER_LOOPS)
    
    mywin.close()
    dataFile.close()
    core.quit()