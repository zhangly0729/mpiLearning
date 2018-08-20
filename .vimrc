"基本配置
set nu "行号
filetype on "文件类型
filetype plugin on
filetype plugin indent on

"Tab 字符
set autoindent 
set tabstop=2
set softtabstop=2
set shiftwidth=2
set expandtab

"VIM配色
syntax enable
set background=dark
"colorscheme solarized

"主题
colorscheme desert

"Ctags
"set tags=/home/echo/src-master/api/tags,/home/echo/src-master/user/tags,./tags;
set tags=./tags;

"折叠
set foldmethod=indent

"SConstruct 高亮
au BufRead,BufNewFile SConstruct set filetype=python
syntax on
syntax enable

"中文显示
set ambiwidth=double
set fileencodings=utf-8,gbk,utf-16le,cp1252,iso-8859-15,ucs-bom
set termencoding=utf-8
set encoding=utf-8

"设备无备份文件
set writebackup
set nobackup

"缩进线条设置
set list
set listchars=tab:\|\ 
"设定文件浏览器目录为当前目录
set autochdir

"设置不自动换行
set nowrap

"显示下方的横向滚动条
set guioptions+=b

"Fortran语言制表符设置
let fortran_have_tabs=1


set tabstop=6
set cindent shiftwidth=6
set autoindent shiftwidth=6

"允许Fortran代码折叠
let fortran_fold=1

"设置代码折叠的方式，这样每个program、module、subroutine、function都可以折叠了
set foldmethod=syntax

"如果前面启用了代码折叠，那么文件一打开代码全部是折叠的，需再按"zO"打开全部折叠的代码
"如果想在文件打开后所有折叠都自动展开，请加入以下配置
set foldlevelstart=99

"设置代码折叠符号（行号左侧），可要可不要，看自己喜欢了
set foldcolumn=4

"设置完毕后重启VIM，打开代码文件，在命令模式下，可以为program、module、subroutine、function折叠代码，常用命令如下：

"zc：折叠代码
"zo：展开代码
"zC：折叠所有代码
"zO：展开所有代码

"map <c-h>,c<space> 


"块代码注释
",ca，在可选的注释方式之间切换，比如C/C++ 的块注释/* */和行注释//
",cc，注释当前行
",c，切换注释/非注释状态
",cs，以”性感”的方式注释
",cA，在当前行尾添加注释符，并进入Insert模式
",cu，取消注释
"Normal模式下，几乎所有命令前面都可以指定行数。  比如  输入  6,cs
"的意思就是以性感方式注释光标所在行开始6行代码
"Visual模式下执行命令，会对选中的特定区块进行注释/反注释
let mapleader=","
set nocompatible                " be iMproved
" Add spaces after comment delimiters by default
 let g:NERDSpaceDelims = 1

"配置WinManager
let g:winManagerWindowLayout='FileExplorer|TagList|BufExplorer'
let g:winManagerWidth=35       "这里设置左侧栏目的宽度
nmap <F3> :WMToggle<cr>      "映射F3键为开关Winmanager

"开启VIM后，自动使用NeoComplete
let g:neocomplcache_enable_at_startup = 1 
"let g:NeoComplCache_EnableAtStartup = 1
