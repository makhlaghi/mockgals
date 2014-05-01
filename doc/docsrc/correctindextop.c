#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int 
main(int argc, char *argv[])
{
  ssize_t read;
  size_t len=0;
  char *line=NULL;
  FILE *in, *out;
  char tmpname[]="purejunk.tmp";

  char initline[]="Next: <a href=\"About.html#About\" accesskey=\"n\" rel=\"next\">About</a>, Previous: <a href=\"../dir/index.html\" accesskey=\"p\" rel=\"prev\">(dir)</a>, Up: <a href=\"../dir/index.html\" accesskey=\"u\" rel=\"up\">(dir)</a> &nbsp; [<a href=\"#SEC_Contents\" title=\"Table of contents\" rel=\"contents\">Contents</a>]</p>\n";

  char changetoline[]="Next: <a href=\"About.html#About\" accesskey=\"n\" rel=\"next\">About</a>, Previous: <a href=\"../index.html\" accesskey=\"p\" rel=\"prev\">(MockGals Home)</a>, Up: <a href=\"../index.html\" accesskey=\"u\" rel=\"up\">(MockGals Home)</a> &nbsp; [<a href=\"#SEC_Contents\" title=\"Table of contents\" rel=\"contents\">Contents</a>]</p>\n";

  in=fopen(argv[1], "r");
  out=fopen(tmpname, "w");
  if (in == NULL) exit(EXIT_FAILURE);
  if (out == NULL) exit(EXIT_FAILURE);

  /* Read until you reach the top line to change */
  while ((read = getline(&line, &len, in)) != -1)
    {
      if(strcmp(line, initline)==0) break;
      fprintf(out, "%s", line);
    }

  /* Put in the correct line. */
  fprintf(out, "%s", changetoline);

  /* Read until you reach the bottom line to change */
  while ((read = getline(&line, &len, in)) != -1)
    {
      if(strcmp(line, initline)==0) break;
      fprintf(out, "%s", line);
    }

  /* Put in the correct line. */
  fprintf(out, "%s", changetoline);

  /* Continue on with the rest of the file: */
  while ((read = getline(&line, &len, in)) != -1)
    fprintf(out, "%s", line);

  /* Free the space allocated for line in getline. */
  if(line) free(line);

  fclose(in);
  fclose(out);

  rename(tmpname, argv[1]);
}
