library(DiagrammeR)

# graphical ---------------------------------------------------------------
#1st graph
grViz(
  "
  digraph dot {
    graph [compound = true, nodesep = .5, ranksep = .25,
           color = black]
      node [shape = circle,style = filled,fillcolor = white,color = black,label = '&pi;'] pi
      node [shape = circle,style = filled,fillcolor = white,color = black,label = 'z@_{1}'] z1
      node [shape = circle,style = filled,fillcolor = white,color = black,label = 'z@_{2}'] z2
      node [shape = circle,style = filled,fillcolor = white,color = white,label = '・・・', fixedsize = True,width=0.3,height=0.3] zk
      node [shape = circle,style = filled,fillcolor = white,color = black,label = 'z@_{n}',fixedsize = False] zn
      edge [color = black]
        pi -> z1
        pi -> z2
        pi -> zk [color = white,style = dotted, arrowhead = none,penwidth = .1]
        pi -> zn
        {rank=min;pi}
        {rank=same;z1;z2;zk;zn}
    }",
  engine = "dot"
)

# graphical ---------------------------------------------------------------
#2nd graph
grViz(
  "
  digraph dot {
    graph [compound = true, nodesep = .5, ranksep = .25,
           color = black]
      node [shape = circle,style = filled,fillcolor = white,color = black,label = '&pi;'] pi
       subgraph cluster1 {
          labelloc=b
          labeljust =r
          label = 'n'
            node [shape = circle,style = filled,fillcolor = white,color = black,label = 'z@_{i}'] zi
        }
      edge [color = black]
        pi -> zi
        {rank=min;pi}
    }",
  engine = "dot"
)

# graphical ---------------------------------------------------------------
#3rd graph
grViz(
  "
  digraph dot {
    graph [compound = true, nodesep = .5, ranksep = .25,
           color = black,newrank=true]
      node [shape = circle,style = filled,fillcolor = white,color = black,label = '&pi;'] pi
      node [shape = circle,style = filled,fillcolor = white,color = black,label = '&eta;'] eta
       subgraph cluster1 {
          labelloc= b
          labeljust = r
          label = 'n'
            node [shape = circle,style = filled,fillcolor = white,color = black,label = 'z@_{i}'] zi
            node [shape = circle,style = filled,fillcolor = white,color = black,label = 'x@_{i}'] xi
        edge [color = black]
        zi -> xi
       }
        subgraph cluster2{
          labelloc = b
          labeljust = r
          label = 'K'
          node [shape = circle,style = filled,fillcolor = white,color = black,label = '&phi;@_{i}'] phii
        }
      edge [color = black]
        pi -> zi
        phii -> xi
        eta -> phii
  {rank=min;pi}
  {rank=same;xi,phii}
    }",
  engine = "dot"
)
# graphical(LDA) ---------------------------------------------------------------
#
grViz(
  "
  digraph dot {
    graph [compound = true, nodesep = .5, ranksep = .25,
           color = black,newrank=true,rankdir=LR]
      node [shape = circle,style = filled,fillcolor = white,color = black,label = '&alpha;'] alpha
      node [shape = circle,style = filled,fillcolor = white,color = black,label = '&beta;'] beta
       subgraph clusterM {
          labelloc= b
          labeljust = r
          label = 'M'
          node [shape = circle,style = filled,fillcolor = white,color = black,label = '&theta;@_{d}'] thetad
          edge [color = black]
          subgraph clustern{
            labelloc= b
            labeljust = r
            label = 'n@_{d}'
            node [shape = circle,style = filled,fillcolor = white,color = black,label = 'z@_{d,i}'] zdi
            node [shape = circle,style = filled,fillcolor = white,color = black,label = 'w@_{d,i}'] wdi
            edge [color = black]
          }  
        }
        subgraph clusterK{
          labelloc = b
          labeljust = r
          label = 'K'
          node [shape = circle,style = filled,fillcolor = white,color = black,label = '&phi;@_{k}'] phik
        }
      edge [color = black]
        alpha -> thetad
        beta -> phik
        thetad -> zdi
        zdi -> wdi
        phik -> wdi
        {rank=min;beta}
        {rank=same;wdi,phik}
    }",
  engine = "dot"
)




