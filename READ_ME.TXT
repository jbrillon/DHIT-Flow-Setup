This program (source text box.for and file box.inc
included into source text through INCLUDE directive)
can be used both to generate a 3D velocity field
with a prescribed energy spectrum, E(k), and to
extract an energy spectrum from a given 3D velocity
field.

To run the program:
1) produce an executable box.exe
2) run it by typing

             box.exe

The program will read the input information from
the file with a fixed name box.in.

A sample file box.in which we are sending contains 
sufficient comments at the right part of each line 
and also at the end of the file.

In addition, we are sending the following files:

   - file e_exp.tec (in Tecplot format) with the
     non-dimensional experimental energy spectra 
     E(k) from the paper by Comte-Bellot and Corsin 
     at time-moments t=0, t=0.87, and t=2.0 
     (3 Tecplot zones for 3 moments).

   - file energy.prf with initial (at t=0) energy
     spectrum from the paper of Comte-Bellot and
     Corsin E(k), k=0, 1, 2,... (the data for those
     k-numbers is obtained by interpolation of the
     experimental data from the 1st zone of the file
     e_exp.tec). This file is an input file for the
     program in the "generate velocity field" mode.

   - file spec.tec with a sample output spectrum
     obtained as follows:
        - running the program in the "generate velocity
          field" mode with the file energy.prf as an input
          and unformatted file velocity.fld as an
          output;
        - switching to "extract energy spectrum" mode
          and running the program with the generated
          at the previous stage file velocity.fld as
          an input and file spec.tec as an output.
          Note that E(k) values in the file spec.tec
          coincide with the first E(k) values of the
          file energy.prf.
