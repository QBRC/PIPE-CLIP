#!/usr/bin/python
from Tkinter import *
import tkFileDialog
import threading
import re
from tkFileDialog import askopenfilename
import tkMessageBox

root = Tk()
root.title("PIPECLIP GUI")

bamFileEntryString = StringVar()
prefixEntryString = StringVar()
matchLenString = StringVar()
maxMismatchString= StringVar()
fdrEnrichedClusterString = StringVar()
fdrReliableMutationString = StringVar()

def BAMFileSelectionCallback():
  file = askopenfilename(filetypes = [("Bam Files","*.BAM")])
  bamFileEntryString.set(file)

def PreBAMFileSelectionCallback():
  thread = threading.Thread(target=BAMFileSelectionCallback)
  thread.start()

def popupError(msg):
  tkMessageBox.showinfo("Parameter Error",msg)

row = 0
col = 0


prefixLabel = Label(root,text="Output Prefix").grid(row=row,column=0)
prefixEntry = Entry(root, bd=3, textvariable=prefixEntryString).grid(row=row,column=1)
row += 1

matchLenLabel = Label(root,text="Shorted matched segment length").grid(row=row,column=0)
matchLenEntry = Entry(root,bd=3,textvariable=matchLenString).grid(row=row,column=1)
row += 1

maxMismatchLabel = Label(root,text="Maximum mismatch number").grid(row=row,column=0)
maxMismatchEntry = Entry(root,bd=3,textvariable=maxMismatchString).grid(row=row,column=1)
row += 1

pcrRemovalLabel = Label(root,bd=3,text="PCR Removal").grid(row=row,column=0)
pcrRemovalListbox = Listbox(root,selectmode="single",height=3,exportselection=0)
pcrRemovalListbox.insert(0,"No removal")
pcrRemovalListbox.insert(1,"Same-start removal")
pcrRemovalListbox.insert(2,"Same-seq removal")
pcrRemovalListbox.grid(row=row,column=1)
row += 1

fdrEnrichedLabel = Label(root,text="FDR for enriched clusters").grid(row=row,column=0)
fdrEnrichedEntry = Entry(root,bd=3,textvariable=fdrEnrichedClusterString).grid(row=row,column=1)
row += 1

clipTypeLabel = Label(root,bd=3,text="CLIP Type").grid(row=row,column=0)
clipTypeListbox = Listbox(root,selectmode="single",height=4,exportselection=0)
clipTypeListbox.insert(0,"HITS-CLIP")
clipTypeListbox.insert(1,"PAR-4SU")
clipTypeListbox.insert(2,"PAR-6SG")
clipTypeListbox.insert(3,"iClip")
clipTypeListbox.grid(row=row,column=1)
row += 1

fdrReliableMutationLabel= Label(root,text="FDR for reliable mutations").grid(row=row,column=0)
fdrReliableMutationEntry = Entry(root,bd=3,textvariable=fdrReliableMutationString).grid(row=row,column=1)
row += 1

speciesLabel = Label(root,bd=3,text="Species").grid(row=row,column=0)
speciesListbox = Listbox(root,selectmode="single",height=3,exportselection=0)
speciesListbox.insert(0,"mm10")
speciesListbox.insert(1,"hg19")
speciesListbox.insert(1,"mm9")
speciesListbox.grid(row=row,column=1)
row += 1

bamFileButton = Button(root,text="Select BAM File",command = PreBAMFileSelectionCallback).grid(row=row,column=0)
bamFileEntry = Entry(root, bd=3,textvariable=bamFileEntryString,state=DISABLED).grid(row=(row),column=1)
row += 1

def processCommandArgs():
  infile =  bamFileEntryString.get()
  outputPrefix = prefixEntryString.get()
  matchLength = matchLenString.get()
  mismatch = maxMismatchString.get()
  fdrCluster = fdrEnrichedClusterString.get()
  fdrMutation = fdrReliableMutationString.get()
  error = False
  if not len(pcrRemovalListbox.curselection()) > 0:
    if not error:
      popupError("Please chhose a PCR Removal")
    error = True
  else:
    pcr =  pcrRemovalListbox.curselection()[0]

  if not len(clipTypeListbox.curselection()) > 0:
    if not error:
      popupError("Please choose a clip type")
    error = True
  else:
    clipType = clipTypeListbox.curselection()[0]

  if not len(speciesListbox.curselection()) > 0:
    if not error:
      popupError("PLease selection a species")
    error = True
  else:
    species =  speciesListbox.curselection()[0]

  if not re.match('[1-4]',mismatch,flags=0):
    if not error:
      popupError("Invalid mismatch parameter; Value Range: [1,4]")
    error = True

  if not re.match('\d+',matchLength,flags=0):
    if not error:
        popupError("Invalid match length parameter")
    error = True
  else:
    if(int(matchLength) < 10):
      popupError("Expected range of values for match length parameter: [10,)");
      error = True

  if not re.match('0?\.\d+',fdrCluster,flags=0):
    if not error:
      popupError("Invalid fdr cluster; Expected Value: (0,1)")
    error = True

  if not re.match('0?\.\d+',fdrMutation,flags=0):
    if not error:
      popupError("Invalid fdr mutation parameter; Expected Value: (0,1)")
    error = True

  if not re.search('\.bam$',infile,flags=0):
    if not error:
      popupError("Invalid BAM File")
    error = True

  if not re.match('^[a-zA-Z0-9\/]+$',outputPrefix,flags=0):
    if not error:
      popupError("Invalid output prefix, must follow regex: ^[a-zA-Z0-9\/]+$")
    error = True

  if not error:
    print "Running analysis"
    import pipeclip
    pipeclip.runPipeClip(infile,outputPrefix,matchLength,mismatch,pcr,fdrCluster,clipType,fdrMutation,species)
    print "Done with analysis"

def runGui(runButton):
  runButton.configure(state=NORMAL,text="Run")

runButton = Button(root,text="Loading PipeClip library...",state=DISABLED,command=processCommandArgs)
runButton.grid(row=row,column=0)
t = threading.Thread(target=runGui,args=(runButton,))
t.start()
root.mainloop()
t.join()
