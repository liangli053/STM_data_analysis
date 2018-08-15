Scanning tunneling microscopy (STM) + Atomic force microscopy (AFM) recorded the changes in the <br/>
frequency of the tip when scanning over a molecule at a certain height Z (in angstrom). <br/>
This code extract the useful part of the data (upper plateaus), and then smooth the pleateaus. <br/>
Integrate the frequency delta_f over Z twice to get the potential. Then compute the derivative of <br/>
the potential with regard to scanning range to get the lateral force required to move the molecule. <br/>
Raw experimental data were stored in the folder raw_data. <br/>
To compute force, run ./Force_2move_atom.py
