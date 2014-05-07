The docs:
=========

The sources for making the documentation are in `./docsrc`.

To compile the documentation you need to run `run.sh`. But `run.sh' 
depends on two compiled files (from `docsrc/`), to compile them (only
once is enough), simply run the `makerunrequirements.sh` that is 
available in this directory. It will compile these two programs that
will make two changes to the HTML version of this manual.

* `addscripttohtmlhead`: This program will put a script in the header
   of the HTML file (currently for MathJax to view equations and
   Google Analytics to analyse traffic). To make this program, you can
   simply run this command in this 'doc' directory: 
   `$ gcc -o addscripttohtmlhead ./docsrc/addscripttohtmlhead.c`.

* `correctindextop`: By default, Texinfo will put a `(dir)` link on
   the top and bottom of the `index.html` that it produces.  The
   address it provides does not fit my webpage! So This program will
   change them to [MockGal's
   Webpage](http://astr.tohoku.ac.jp/~akhlaghi/mockgals.html).  You
   can compile it with 
   `$ gcc -o correctindextop ./docsrc/correctindextop.c`

Leave the two compiled programs in this directory, they will be 
needed and don't need frequent change!
