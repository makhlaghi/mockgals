#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void
addtexttohtml(char *scriptname, char *htmlname)
{
  ssize_t read;
  size_t len=0;
  char *line=NULL;
  FILE *in, *out, *script;
  char tmpname[]="purejunk.tmp";

  in=fopen(htmlname, "r");
  out=fopen(tmpname, "w");
  script=fopen(scriptname, "r");
  if (in == NULL) exit(EXIT_FAILURE);
  if (out == NULL) exit(EXIT_FAILURE);
  if (script == NULL) exit(EXIT_FAILURE);

  /* Read until you reach the </head> tag. */
  while ((read = getline(&line, &len, in)) != -1)
    {
      if(strcmp(line, "</body>\n")==0) break;
      fprintf(out, "%s", line);
    }

  /* Put in the script from the script file. */
  while ((read = getline(&line, &len, script)) != -1)
    fprintf(out, "%s", line);
  fprintf(out, "\n\n\n</body>\n\n\n");

  /* Continue on with the rest of the file: */
  while ((read = getline(&line, &len, in)) != -1)
    fprintf(out, "%s", line);

  /* Free the space allocated for line in getline. */
  if(line)
    free(line);

  fclose(in);
  fclose(out);
  fclose(script);

  rename(tmpname, htmlname);
}


/* Call this program like this: 
a.out htmlheaderscript.txt ./HTMLFOLDER/*.html  */
int
main(int argc, char *argv[])
{
  int i;

  for(i=2;i<argc;i++)
    addtexttohtml(argv[1], argv[i]);

  return 0;
}
