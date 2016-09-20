%
%   Load expl3 package.
\RequirePackage{expl3}[2013/07/28]
\@ifpackagelater{expl3}{2013/07/28} {} {%
    \PackageError{UWMadThesis}{Version of l3kernel is too old}
      {%
        Please install an up to date version of l3kernel\MessageBreak
        using your TeX package manager or from CTAN.
      }%
    \endinput
}%
%
%   Turn on expl3.
\ExplSyntaxOn
%
%   Declare the LaTeX format needed
\NeedsTeXFormat{LaTeX2e}
%
%   Provide the class
\ProvidesClass{UWMadBeamer}
    [2014/05/10~v0.1~UW-Madison~presentation~class]
%
%   Pass all user options directly to beamer
\DeclareOption*{
    \PassOptionsToClass{\CurrentOption}{beamer}
}
%
%   Declare defaults and process the options
\ExecuteOptions{12pt,t,xcolor=table}
\ProcessOptions
%
%   load beamer
\LoadClass{beamer}[2013/01/01]
%
%   Load packages
\RequirePackage{xparse}
\RequirePackage{appendixnumberbeamer}
\RequirePackage[T1]{fontenc}
\RequirePackage{lmodern}
\RequirePackage{microtype}
\RequirePackage{graphicx}
\RequirePackage{xparse}
%
%   Define the layout and color themes
\useoutertheme[subsection=false]{miniframes}
\useinnertheme{circles}
\usecolortheme{UWMad}
\hypersetup{colorlinks,linkcolor=}
\let\HrefOld\href
\DeclareDocumentCommand \href { m m } {
    \textcolor{UWMadRed}{\HrefOld{#1}{#2}}
}
\hypersetup{urlcolor=UWMadRed}
%\useoutertheme{default}
%\useinnertheme{circles}
%
%   Modify the color scheme a bit
%
%   Remove navigation symbols
\setbeamertemplate{navigation~symbols}{}
%
\setbeamertemplate{headline}{
    \begin{beamercolorbox}[wd=\paperwidth]{section~in~head/foot}
        \vskip2pt\insertnavigation{\paperwidth}\vskip2pt
    \end{beamercolorbox}%
    \begin{beamercolorbox}[wd=\paperwidth,ht=0.004\paperheight]{fine~separation~line}
    \end{beamercolorbox}%
}
%
%
%
%   Footline: Date - Title - Slide Count
\setbeamertemplate{footline}{
    %
    %   Vertical corrections and gold box separator
    \begin{beamercolorbox}[wd=\paperwidth,ht=0.004\paperheight]{fine~separation~line}
    \end{beamercolorbox}%
    %
    %   Encompassing box
    \begin{beamercolorbox}{palette~tertiary}
        %
        %   Left, small box for date
        \begin{beamercolorbox}
            [wd=0.1\paperwidth,ht=0.024\paperheight,
            dp=0.00266666\paperheight,leftskip=0.01\paperwidth]
            {palette~tertiary}
            \insertshortdate{}
        \end{beamercolorbox}
        %
        %   Center, large box for title
        \begin{beamercolorbox}
            [center,wd=0.7999\paperwidth,
            ht=0.024\paperheight,dp=0.00266666\paperheight]
            {palette~tertiary}
            \insertshorttitle{}
        \end{beamercolorbox}
        %
        %   Right, small box for slide count
        \begin{beamercolorbox}
            [right,wd=0.1\paperwidth,ht=0.024\paperheight,
            dp=0.00266666\paperheight,rightskip=0.01\paperwidth]
            {palette~tertiary}
            \insertframenumber/\inserttotalframenumber
        \end{beamercolorbox}
        %
        %   Bit ox extra vertical space
        \vbox{}\vskip0.01\paperheight
    \end{beamercolorbox}
}
%
%
%
%
%
%
%
%
\tl_new:N \g_UWMadBeamer_Title_Short_tl
\tl_new:N \g_UWMadBeamer_Title_Long_tl
\tl_new:N \g_UWMadBeamer_Author_Short_tl
\tl_new:N \g_UWMadBeamer_Author_Long_tl
\tl_new:N \g_UWMadBeamer_Institution_Short_tl
\tl_new:N \g_UWMadBeamer_Institution_Long_tl
\tl_new:N \g_UWMadBeamer_Department_Short_tl
\tl_new:N \g_UWMadBeamer_Department_Long_tl
\tl_new:N \g_UWMadBeamer_Program_Short_tl
\tl_new:N \g_UWMadBeamer_Program_Long_tl
\tl_new:N \g_UWMadBeamer_Logo_FileName_tl
%
%
\cs_new_eq:NN \__UWMadBeamer_Author:nn     \author
\cs_new_eq:NN \__UWMadBeamer_Title:nn      \title
\cs_new_eq:NN \__UWMadBeamer_Institute:nn  \institute
\cs_undefine:N \author
\cs_undefine:N \title
\cs_undefine:N \institute
%
%
% Identification commands
\DeclareDocumentCommand \author { o m } {
    \IfNoValueTF {#1} {
        \tl_gset:Nn \g_UWMadBeamer_Author_Short_tl {#2}
        \tl_gset:Nn \g_UWMadBeamer_Author_Long_tl  {#2}
        \__UWMadBeamer_Author:nn[#2]{#2}
    } {
        \tl_gset:Nn \g_UWMadBeamer_Author_Short_tl {#1}
        \tl_gset:Nn \g_UWMadBeamer_Author_Long_tl  {#2}
        \__UWMadBeamer_Author:nn[#1]{#2}
    }
}
\DeclareDocumentCommand \institute { o m } {
    \IfNoValueTF {#1} {
        \tl_gset:Nn \g_UWMadBeamer_Institute_Short_tl {#2}
        \tl_gset:Nn \g_UWMadBeamer_Institute_Long_tl  {#2}
        \__UWMadBeamer_Institute:nn[#2]{#2}
    } {
        \tl_gset:Nn \g_UWMadBeamer_Institute_Short_tl {#1}
        \tl_gset:Nn \g_UWMadBeamer_Institute_Long_tl  {#2}
        \__UWMadBeamer_Institute:nn[#1]{#2}
    }
}
\DeclareDocumentCommand \title { o m } {
    \IfNoValueTF {#1} {
        \tl_gset:Nn \g_UWMadBeamer_Title_Short_tl {#2}
        \tl_gset:Nn \g_UWMadBeamer_Title_Long_tl  {#2}
        \__UWMadBeamer_Title:nn[#2]{#2}
    } {
        \tl_gset:Nn \g_UWMadBeamer_Title_Short_tl {#1}
        \tl_gset:Nn \g_UWMadBeamer_Title_Long_tl  {#2}
        \__UWMadBeamer_Title:nn[#1]{#2}
    }
}
\DeclareDocumentCommand \department { o m } {
    \IfNoValueTF {#1} {
        \tl_gset:Nn \g_UWMadBeamer_Department_Short_tl {#2}
        \tl_gset:Nn \g_UWMadBeamer_Department_Long_tl  {#2}
    } {
        \tl_gset:Nn \g_UWMadBeamer_Department_Short_tl {#1}
        \tl_gset:Nn \g_UWMadBeamer_Department_Long_tl  {#2}
    }
}
%
%
%   Display outline at start of every (sub)section.
\AtBeginSection[]{
    \begin{frame}<beamer>
        \frametitle{Outline}
        \tableofcontents[currentsection]
    \end{frame}
}
\AtBeginSubsection[]{
    \ifnum\value{subsection}>1
    \begin{frame}<beamer>
        \frametitle{Outline}
        \tableofcontents[currentsection,currentsubsection]
    \end{frame}
    \fi
}
\setlength{\parskip}{0.5em}
%
%
%
%
\ExplSyntaxOff