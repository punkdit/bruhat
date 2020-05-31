#!/usr/bin/env python3

import os
import hashlib
from subprocess import Popen, PIPE

from bruhat.render import load_svg

def file_exists(name):
    try:
        os.stat(name)
        return True
    except:
        pass
    return False


def verbose_command(cmd):
    ret = os.system(cmd)
    assert ret == 0, "%r failed with return value %d"%(cmd, ret)


def command(cmd):
    p = Popen(cmd, shell=True,
              stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
    (child_stdin,
     child_stdout,
     child_stderr) = (p.stdin, p.stdout, p.stderr)
    ret = p.wait()
    assert ret == 0, "%r failed with return value %d"%(cmd, ret)

#command = verbose_command


def tex_output(text):
    lines = []
    lines.append(r"\def\folio{}")
    lines.append(r"%s\bye"%text)
    return '\n'.join(lines)


def latex_output(text):
    lines = []
    lines.append(r"\documentclass[11pt]{article}")
    lines.append(r"\pagenumbering{gobble}")
    lines.append(r"\usepackage{fontspec}")
    lines.append(r"\setromanfont{Gentium Book Basic}")
    lines.append(r"\def\folio{}")
    lines.append(r"\begin{document}")
    lines.append(text)
    lines.append(r"\end{document}")
    return '\n'.join(lines)


def make_text(text, tex_engine="pdftex"):
    assert tex_engine in "pdftex xetex xelatex pdflatex".split()
    cache = "__bruhat__"
    if not file_exists(cache):
        os.mkdir(cache)
    os.chdir(cache) # <---------- chdir <-----

    if tex_engine == "pdftex" or tex_engine == "xetex":
        output = tex_output(text)
    elif tex_engine == "xelatex" or tex_engine=="pdflatex":
        output = latex_output(text)
    else:
        assert 0, tex_engine

    data = output.encode('utf-8')
    stem = hashlib.sha1(data).hexdigest()

    tex_name = "%s.tex"%stem
    svg_name = "%s.svg"%stem
    pdf_name = "%s.pdf"%stem

    if not file_exists(svg_name):

        f = open(tex_name, 'w')
        print(output, file=f)
        f.close()
    
        command("%s %s"%(tex_engine, tex_name))
        command("pdf2svg %s %s" % (pdf_name, svg_name))

    item = load_svg.load(svg_name)

    os.chdir("..") # <---------- chdir <-----

    return item


def test():

    item = make_text("hi there!")
    print(item)


if __name__ == "__main__":

    test()

