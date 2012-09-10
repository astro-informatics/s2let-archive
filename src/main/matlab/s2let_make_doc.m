function s2let_make_doc(s2letpath)

cd(s2letpath)
m2html('mfiles', 'src/main/matlab', 'htmldir', 'doc/matlab');

end