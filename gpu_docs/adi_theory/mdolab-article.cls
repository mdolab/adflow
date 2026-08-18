%%
%% This is file `mdolab-article.cls'
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{mdolab-article}

% Define an option called nopub
\newif\if@published
\@publishedtrue % default is published
\DeclareOption{nopub}{\@publishedfalse}

% Define an option called noinfobox
\newif\if@infobox
\@infoboxtrue % default is to include the info box
\DeclareOption{noinfobox}{\@infoboxfalse}

% Define an option called preprint
\newif\if@preprint
\@preprintfalse % default is postprint
\DeclareOption{preprint}{\@preprinttrue}

% Define an option called defaultfonts that will use the article class fonts
% This also allows the user to set their own font without causing package conflicts
\newif\if@mdobookfonts
\@mdobookfontstrue % default is to use mdobook fonts
\DeclareOption{defaultfonts}{\@mdobookfontsfalse}

% Pass the rest of class options to the underlying article class
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions

% Base it of the default article class
\LoadClass[11pt]{article}

% Redefine the page margins
\RequirePackage[left=1in,right=1in,top=1in,bottom=1in]{geometry}

\RequirePackage[
    pdftex,
    bookmarks=true,
    bookmarksnumbered=true,
    colorlinks=true,
    linkcolor=qblue,
    citecolor=qblue,
    urlcolor=qblue,
]{hyperref}

% To modify the heading styles more thoroughly use the titlesec package
\RequirePackage{titlesec}

% Other packages to load by default
\RequirePackage{amsmath}
\RequirePackage{graphicx}
\RequirePackage[bf,footnotesize]{caption}
\RequirePackage{indentfirst}
\RequirePackage{afterpage}

% Set the mdobook fonts, some settings from kaobook
\if@mdobookfonts
\RequirePackage[scaled=.97,helvratio=.93,p,theoremfont]{newpxtext} % Palatino for serif
\RequirePackage[vvarbb,smallerops,bigdelims]{newpxmath} % Palatino for math
\linespread{1.07} % More space between lines is recommended for Palatino
\RequirePackage[defaultsans]{opensans} % Open Sans for sans serif
\RequirePackage[scaled=0.85]{beramono} % Monospace font
\else
\RequirePackage{amssymb} % This conflicts with newpxmath
\fi

% Author block
\RequirePackage[blocks]{authblk}
\renewcommand{\Authfont}{\normalsize\upshape}
\renewcommand{\Affilfont}{\normalsize\itshape}
\setlength{\affilsep}{0pt}

% define some colors
\RequirePackage[table,xcdraw]{xcolor}
\colorlet{qblue}{blue!50!black}

% The fields to be used for the title page
\newcommand{\displaymdolaburl}{mdolab.engin.umich.edu}
\newcommand{\fullmdolaburl}{https://mdolab.engin.umich.edu}
\newcommand{\@officialurl}{Insert preprint URL!}
\newcommand{\@papercitation}{Insert paper citation!}
\newcommand{\@keywords}{Insert paper keywords!}

% Commands to set the titlepage fields
\newcommand{\officialurl}[1]{\renewcommand{\@officialurl}{#1}}
\newcommand{\papercitation}[1]{\renewcommand{\@papercitation}{#1}}
\newcommand{\keywords}[1]{\renewcommand{\@keywords}{#1}}

% Commands for the user to be able to use the defined fields.
\newcommand{\insertofficialurl}{\@officialurl}
\newcommand{\insertpapercitation}{\@papercitation}
\newcommand{\insertkeywords}{\@keywords}

% Command to say postprint or preprint
\if@preprint
\newcommand{\printtype}{preprint}
\else
\newcommand{\printtype}{postprint}
\fi

\definecolor{mdoDarkBlue}{HTML}{003268}

% Command to make stylized ribbon for preprint/postprint description
\RequirePackage{tcolorbox}
\newtcolorbox{headerColorBox}{
    colback=mdoDarkBlue,
    colframe=mdoDarkBlue,
    coltext=white,
    grow to left by=1in,
    grow to right by=0.3in,
    left*=0mm,
    sharp corners,
}

% Blue box with paper information
\newcommand{\makeinfobox}{
  \begin{center}
      % Published paper
      \begin{headerColorBox}
      {\bfseries \small
          This is a \printtype{} of the following article:} \\
      \noindent
      {\small\insertpapercitation} \\
      {\bfseries \small
      \if@published
      The published article may differ from this \printtype{} and is available at: \\
      {\hypersetup{urlcolor=white}\url{\insertofficialurl}.}
      \else
      The published article may differ from this \printtype{}.
      \fi
      }
      \end{headerColorBox}
  \end{center}
}

% Adjust the title design
\newcommand{\maketitlepage}{%
  \begingroup
  \centering
  {\LARGE\sffamily\bfseries \@title \par}
  \vskip\baselineskip
  \@author\par
  \endgroup
  \vskip2\baselineskip
}
% This some before-and-after code that surrounds the title page.  It shouldn't need to be modified.
% I've pulled out the part the actually typesets the title page and placed it in the \maketitlepage command above.
\renewcommand\maketitle{
  \if@infobox
    \newgeometry{left=1in,right=1in,top=0.5in,bottom=1in}
    \afterpage{\aftergroup\restoregeometry}
    \makeinfobox
    \vspace{0.1in}
  \fi
  \maketitlepage%
  \setcounter{footnote}{0}
  % update hyperref with PDF metadata since we have that info right here
  \hypersetup{
    % note that the macro \AB@author is taken from the authblk package which redefines \@author to contain the affiliation
    % which cannot be used here, but they offer \AB@author which contains only the authors
    pdfauthor={\AB@author}, % the authors
    pdftitle={\@title}, % paper title
    pdfkeywords={\insertkeywords}, % keywords, like those used in journal papers for subject classification
  }
}

%% Abstract title: Set to Sans Serif font
%\renewenvironment{abstract}{%
%      \@beginparpenalty\@lowpenalty
%      \begin{center}%
%      \vskip 0\p@
%      \bfseries\sffamily \abstractname
%      \@endparpenalty\@M
%      \end{center}
%      \vskip -10\p@
%      }%

% Define `abstract' environment
% \def\abstract{\topsep=0pt\partopsep=0pt\parsep=0pt\itemsep=0pt\relax
% \trivlist\item[\hskip\labelsep
% {\Large\bfseries\sffamily\abstractname}]\if!\abstractname!\hskip-\labelsep\fi}
% \def\abstract{\topsep=0pt\partopsep=0pt\parsep=0pt\itemsep=0pt\relax
% \trivlist\item[\hskip\labelsep
% {\large\bfseries\sffamily\abstractname}]\if!\abstractname!\hskip-\labelsep\fi}
\renewenvironment{abstract}{%
%\titlepage
%\null\vfil
\@beginparpenalty\@lowpenalty
\begin{center}%
  \bfseries \sffamily \Large \abstractname
  \@endparpenalty\@M
\end{center}}%
{\par}
%{\par\vfil\null\endtitlepage}

% Section titles: Set to Sans Serif font

\renewcommand\section{\@startsection {section}{1}{\z@}%
    {-3.5ex \@plus -0.5ex \@minus -.2ex}%
    {1.2ex \@plus.1ex}%
    {\sffamily\Large\bfseries}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
    {-3.25ex\@plus -0.5ex \@minus -.2ex}%
    {0.8ex \@plus .2ex}%
    {\sffamily\large\bfseries}}
\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
    {-3.25ex\@plus -0.5ex \@minus -.2ex}%
    {0.8ex \@plus .2ex}%
    {\sffamily\normalsize\bfseries}}

% JM: Below are the original preprint settings, modified from the standard article.cls. I changed them to more reasonable spacing to improve readability (see above)
% \renewcommand\section{\@startsection {section}{1}{\z@}%
%     {-3.5ex \@plus -0.1ex \@minus -.2ex}%
%     {0.2ex \@plus.01ex}%
%     {\sffamily\large\bfseries}}
% \renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
%     {-3.25ex\@plus -0.1ex \@minus -.2ex}%
%     {0.2ex \@plus .01ex}%
%     {\sffamily\normalsize\bfseries}}
% \renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
%     {-3.25ex\@plus -0.1ex \@minus -.2ex}%
%     {0.2ex \@plus .01ex}%
%     {\sffamily\MakeUppercase\small\bfseries}}

% Table of contents: Set to Sans Serif font and tweak the spacing
\setcounter{tocdepth}{1} % set to 2 to get next level
\RequirePackage{tocloft}
\setlength\cftparskip{0pt}
\setlength\cftbeforesecskip{1pt}
\setlength\cftaftertoctitleskip{2pt}
\renewcommand{\cfttoctitlefont}{\Large\bfseries\sffamily}
\renewcommand{\cftsecfont}{\bfseries\sffamily}
\renewcommand{\cftsecleader}{\bfseries\sffamily\cftdotfill{\cftdotsep}}
\renewcommand{\cftsecpagefont}{\bfseries\sffamily}


% Abstract


\endinput
