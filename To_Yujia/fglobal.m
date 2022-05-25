function o=fglobal(gr,pos)
  o=2.^(0.66.*(1.6-gr./log(2))).*2.^(gr.*(0.78+0.15./gr).*(1-pos)+0.25);
end