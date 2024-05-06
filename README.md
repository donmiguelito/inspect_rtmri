# inspect_rtmri
Explore and analyze real-time MRI data

DESCRIPTION
-----------
INSPECT_rtMRI is a Matlab GUI designed to provide an integrated interface for inspection and analysis of companion audio and video files recorded during dynamic MRI sessions.


REQUIREMENTS
------------
- Matlab v.2007a or later
- interface has not been optimized for OS-X, and parts may poorly on Mac versions Matlab
- interface has not been updated to comply with all recent changes in Matlab


QUICK START
-----------
- update Matlab path to include directory </rtMRI_GUI>
- set Matlab Current Folder to directory </rtMRI_GUI>
- launch interface from Matlab command window:
  inspect_rtmri('kaLi')


GENERAL USE
-----------
- the GUI will work with any rtMRI data in *.avi format 
- optimal use requires synchronized companion audio files in *.wav format 
- the GUI is designed to operate in a working directory containing /avi and /wav folders
- the GUI expects to find pairs of audio and video files sharing the same base filename
- initiate the GUI with a single string input argument corresponding to the shared base filename

e.g. dat = inspect_rtmri('file0001');
     will load video data from <./avi/file0001.avi>
           and audio data from <./wav/file0001.wav>


---------------
Michael Proctor
26/06/15
mike.i.proctor@gmail.com

