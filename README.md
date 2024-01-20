# EOS
EOS for GAS and Water

This is a collection of Equation of State for GAS and Water written by others, for a detail explanation refer to Github pages
   - AGA 8 - Written by Eric W. Lemmon
   - GERG 2008  - Written by Eric W. Lemmon
  - Generalized Cubic - J. Res. NIST, 2016
  - IF97 - Ian Bell and the CoolProp team

I collect the codes in a simple application to use as a handy tool ready to use quickly, session is recorder in a log file.

Application start asking name of file to log the current section, extension ".log" is added to name.
Next main menu is shown

<img width="229" alt="image" src="https://github.com/markpicci/EOS/assets/144413593/ba57f05d-8b50-4a9f-9bb1-27b10ccaf8e3">

Selecting "9" for water EOS IAPWS 97 only pressure and temperature are required.

Other EOS for gas requires to define the mixture selecting "m"

<img width="250" alt="image" src="https://github.com/markpicci/EOS/assets/144413593/4c004334-f37f-4df2-ae5b-d4bd024be2b4">

Item g - show available mixture in working directory
Item l - load a saved mixture
Item s - save mixture
Item i - input a new mixture
Item m - show active mixture
Item c - show all mixture that are loaded and select active
Item d - deleet ALL mixture
item e - exit and return to main menu

I would like to note that binary interaction parameter for the Generalized Cubic Peng Robinson equation are quite rude and shall be revised, any contribution is most then welcome.

Licensed under the [GNU Lesser General Public License, Version 3.0 without any warranty or liability.
