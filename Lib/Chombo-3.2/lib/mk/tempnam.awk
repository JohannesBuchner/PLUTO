BEGIN {
  buffer = 0;
}
# gets rid of annoying GNU compiler link-warning 
/`tempnam' is dangerous/ {
  buffer = 0;
  next;
}
/`tmpnam' is dangerous/ {
  buffer = 0;
  next;
}
# filters out annoying linker warnings on AIX machine (seaborg)
/WARNING: Duplicate symbol/ {
  buffer = 0;
  next;
}
#yet another warning that makes no sense
/warning: invalid offsetof from non-POD type/ {
  buffer = 0;
  next;
}
#yet another warning that makes no sense
/pointer to member instead/ {
  buffer = 0;
  next;
}

buffer == 0 {
  save = $0;
  buffer = 1;
  next;
}
buffer == 1 {
  print save;
  save = $0;
}
END {
  if (buffer == 1) {
    print save;
  }
}
