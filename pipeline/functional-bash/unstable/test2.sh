# escapestr_sed()
# read a stream from stdin and escape characters in text that could be interpreted as
# special characters by sed
escape_sed() {
 sed \
  -e 's/\//\\\//g' \
  -e 's/\&/\\\&/g'
}
