'''
Generates images of magnetic field lines and wires 
from data files created using calculate.py

Created on Dec 7, 2011

@author: Eric Leong
with guidance from Professor Wolf
'''
from mayavi import mlab
from mayavi.tools.helper_functions import plot3d
from tvtk.tools import visual
import glob
import os
import pickle

def draw(name, w=[], l=[], mag=False):
    # Load the data into memory.
    f = open(name, 'r')
    wires, fieldlines = pickle.load(f)
    f.close()
    
    # Display each wire.
    wi = 0
    for wire in wires:
        if len(wire) <= 0:
            continue
        
        wx, wy, wz = zip(*wire) # split x, y, z from w
        scalars = [1] * len(wx)
        if wi < len(w):
            w[wi].trait_set(visible=True)
            w[wi].mlab_source.reset(x=wx, y=wy, z=wz, scalars=scalars)
        else:
            # Use a nice copper color.
            w.append(plot3d(wx, wy, wz, scalars, color=(.72, .45, .2), reset_zoom=False))
        wi += 1
    
    # make the extra wires disappear
    for i in range(wi + 1, len(w)):
        w[i].trait_set(visible=False)

    # find vmax, vmin for this set of data
    vmax = -1;
    vmin = -1;
    for line in fieldlines:
        lx, ly, lz, scalars = zip(*line)
        m = max(scalars)
        if m > vmax or vmax == -1:
            vmax = m
        m = min(scalars)
        if m < vmin or vmin == -1:
            vmin = m
    
    # Display each fieldline
    li = 0
    for line in fieldlines:
        lx, ly, lz, scalars = zip(*line)
        
        if li < len(l):
            l[li].trait_set(visible=True)
            l[li].mlab_source.reset(x=lx, y=ly, z=lz, scalars=scalars)
        else:
            if mag:
                # Colored field lines
                l.append(plot3d(lx, ly, lz, scalars, tube_radius=None, reset_zoom=False, vmax=vmax, vmin=vmin))     
            else:
                # White field lines
                scalars = [1] * len(lx)
                l.append(plot3d(lx, ly, lz, scalars, tube_radius=None, reset_zoom=False, color=(1, 1, 1)))
                
        li += 1
    
    # make the extra lines disappear
    for i in range(li + 1, len(l)):
        l[i].trait_set(visible=False)


if __name__ == "__main__":
    path = ''
    
    # Prepare mlab.
    mlab.options.offscreen = True
    fig = mlab.figure(bgcolor=(0, 0, 0), size=(1280, 720))
    visual.set_viewer(fig)
    
    # Iterate through all the files in this folder.
    for name in glob.glob(os.path.join(path, '*.txt')):
        num = int(name[-8:-4])
        
        if num == 1:
            continue
    
        # Disable rendering for faster speed.
        fig.scene.disable_render = True
        draw(name, mag=True)
        fig.scene.disable_render = False
    
        # Set view and save image.
        
        # line
        mlab.view(azimuth=45, elevation=90, distance=6, reset_roll=True, focalpoint=[0, 0, 0])
        mlab.roll(roll=45)
        
        # toroid
        #mlab.view(azimuth=90, elevation=45, distance=20, reset_roll=True, focalpoint=[0, 0, 0])
        
        mlab.savefig(name[:-3] + 'png')
        
        print name
        
        #mlab.show()
        #raw_input("Press Enter")
