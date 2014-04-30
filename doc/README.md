The docs:
=========

The sources for making the documentation are in './docsrc'.

In order to successfully add the scripts to the <head>...</head>
section of the HTML, 'run.sh' depends on the addtohtmlheader.c 
program to be compiled on your system. To do that you can simply
run this command in this 'doc' directory:

    $ gcc -o addscripttohtmlhead ./docsrc/addscripttohtmlhead.c

Leave the compiled program in this directory, it will be needed 
and doesn't need frequent change!
