# Repository.  

Original human 3D atrium model code. (2003-2013).

# Background.  

This is one of the first human 3D atrium model that was developed. It used a geometry provided by
collaborators in Germany, and simulates electrical activity in the cardiac chambers called atria.
In addition to the 3D FD solver, a representative cell to 3D pipeline is also provided to show
how a multi-scale simulation could be performed.

Some of the outputs that emerged from this program were:  
* Kharche et al. (LNCS peer reviewed int. conf. proceedings). https://link.springer.com/chapter/10.1007/978-3-540-72907-5_14  
* Kharche et al. (journal article). https://www.sciencedirect.com/science/article/pii/S0079610708000709
* Kharche et al. (journal article). https://pubmed.ncbi.nlm.nih.gov/17905416/
* Kharche et al. (journal article). https://physoc.onlinelibrary.wiley.com/doi/full/10.1113/jphysiol.2012.229146

This code was (is being) extensively used to support PhD thesis and research by others.

# Dependencies.  

This code is OpenMP based, and has no other dependencies.

# Install.  

* OpenMP. (https://www.openmp.org/)  
* GNU C  
* Data visualization using ParaView  (https://www.kitware.com/open-source/).

# Sources and data description.

The basic program is provided in the sources (02.3DAtrium2003-2012/cprogs_3d_myOriginalHumanAtrium/*.c).
It will simulate approx. 5 seconds worth of data.  

Further, a ready to use multi-scale pipeline from cell to 2D/3D is also provided (see: FIMH2007/FIMHprograms).
Whereas it was used multiple times, it is now a suitable teaching and training instrument.

# Use.

All the code comes with a Makefile, after which the user is expected to set up the environment for number of threads,
and run the binary in that environment.

I have used the above codes in several articles, the first of which is:
Kharche et al. LNCS (2007); 4466: 129–138.  

For others, please see list of publications upto 2012.

# Uptake by other users.

The above codes were used to support multiple PhD theses during 2006-2012 (Uni. Manchester) as well as ensuing
journal articles.

# Maintainer.

Whereas this code is provided as is and works as specificed, the PetSc-Sundials based VCPL package or the BeatBox cardiac simulator package are recommended.
As such, this code is not maintained nor supported but may be used as an educational and training instrument.

# Provenance.

Potential new releases of this code are handled by the lead developer and co-investigator Dr. SR Kharche.

# Acknowledements.

This project was generously funded by the UK British Heart Foundation (2006-2009). project.

# Licence.

BSD 3-Clause License

Copyright (c) 2023, Sanjay R. Kharche, Ph.D.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Dr. Sanjay R. Kharche, Ph.D. (App. Maths).  
January 23rd, 2023.  
Email: Sanjay.Kharche@lhsc.on.ca , skharche@uwo.ca , srk.london@yahoo.com .  
Phone: +15198780685 (mobile).  
Website: https://kcru.lawsonresearch.ca/research/srk/index.html  

