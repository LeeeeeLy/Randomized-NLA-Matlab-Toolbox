# Randomized NLA Matlab Toolbox

## main functions:

| Algorithms: | file name             | from paper # | algorithm # or name listed in the paper                    | simple input tested working | 
| ----------- | --------------------- | ------------ | ---------------------------------------------------------- | --------------------------- | 
| 1           | randRF                |[[1]](#1)     | 4.1                                                        | yes                         |                          
| 2           | Adaptive\_randRF      |[[1]](#1)     | 4.2                                                        | yes                         |                         
| 3           | randPI                |[[1]](#1)     | 4.3                                                        | yes                         |                          
| 4           | randSI                |[[1]](#1)     | 4.4                                                        | yes                         |                          
| 5           | FastRandRF            |[[1]](#1)     | 4.5                                                        | yes                         |                          
| 6           | DirectEigvalueDecopo  |[[1]](#1)     | 5.3                                                        | yes                         |                          
| 7           | EigvalueDecopoRow     |[[1]](#1)     | 5.4                                                        | yes                         |                          
| 8           | EigvalueDecopoNystrom |[[1]](#1)     | 5.5                                                        | yes                         |                          
| 9           | EigvalueDecopoOnePass |[[1]](#1)     | 5.6                                                        | yes                         |                          
| 10          | FixedRank             |[[1]](#1)     | Proto-Algorithm: Solving the Fixed-Rank Problem            | yes                         |                          
| 11          | RandSVD               |[[1]](#1)     | Prototype for Randomized SVD                               | yes                         |                          
| 12          | AERandSVD             |[[2]](#2)     | ACCURACY ENHANCED RANDOMIZED SVD                           | yes                         |                          
| 13          | AEORandSVD            |[[2]](#2)     | ACCURACY ENHANCED RANDOMIZED SVD (WITH ORTHONORMALIZATION) | yes                         |                          
| 14          | BasicRandSVD          |[[2]](#2)     | RSVD                                                       | yes                         |                          
| 15          | SPRandEVDH            |[[2]](#2)     | SINGLE-PASS RANDOMIZED EVD FOR A HERMITIAN MATRIX          | yes                         |                          
| 16          | SPRandSVD             |[[2]](#2)     | SINGLE-PASS RANDOMIZED SVD FOR A GENERAL MATRIX            | yes                         |                          
| 17          | randPowerMethod       |[[3]](#3)     | 4                                                          | yes                         |                          
| 18          | randomizedLanczos     |[[3]](#3)     | 5                                                          | yes                         |                          


## References
<a id="1">[1]</a> 
Halko, Nathan, Per-Gunnar Martinsson, and Joel A. Tropp. 
Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions.
SIAM review 53.2 (2011): 217-288.

<a id="2">[2]</a> 
Martinsson, Per-Gunnar. 
Randomized methods for matrix computations.
The Mathematics of Data 25 (2019): 187-231.

<a id="3">[3]</a> 
Martinsson, Per-Gunnar, and Joel Tropp.
Randomized numerical linear algebra: Foundations & algorithms.
arXiv preprint arXiv:2002.01387 (2020).
The Mathematics of Data 25 (2019): 187-231.

<a id="4">[4]</a> 
Saibaba, Arvind K. 
Randomized subspace iteration: Analysis of canonical angles and unitarily invariant norms.
SIAM Journal on Matrix Analysis and Applications 40.1 (2019): 23-48.


\chapter{Randomized NLA Tool Box for \sc{MATLAB}}

\paragraph{Code Link:}
As a part of this project, we have created a simple randomized numerical linear algebra Matlab Toolbox which can be found on GitHub 
\href{https://github.com/LeeeeeLy/Randomized-NLA-Matlab-Toolbox}{\color{blue} \myul[blue] {https://github.com/LeeeeeLy/Randomized-NLA-Matlab-Toolbox}}
\hspace{1cm}\newline

In the following we will briefly describe the structure of the Tool Box.
\paragraph{Main Functions:}
\begin{itemize}
\item \textbf{FixedRank.m} - This is an implementation of Proto-Algorithm: Solving the Fixed-Rank Problem from paper \cite{4}.
		
		
		This function takes an input matrix A, a target rank k, and a oversampling parameter p and output a $m\times (k+p)$ matrix Q whose columns are orthonormal and whose range approximates the range of A.
         

\item \textbf{RandSVD.m} - This is implementation of Prototype for Randomized SVD from paper \cite{4}.
		
		
		This function takes an input matrix A, a target rank k, and an exponent q (q = 0,1,2) and approximate a rank-2k factorization $A \approx U\Sigma V^{*}$, where U and V are orthonormal, and $\Sigma$ is nonnegative and diagonal.

\item \textbf{randRF.m} - This is implementation of Algorithm 4.1 from paper \cite{4}.
		
		
This function takes an $m\times n$ matrix A and an integer l, it computes an $m \times l$ orthonormal matrix Q whose range approximates the range of A.
\item \textbf{Adaptive\_randRF.m} - This is implementation of Algorithm 4.2 from paper \cite{4}.
		
		
		This function takes an $m \times n$ matrix A, a tolerance $\epsilon$, and an integer r (e.g., r = 10), it computes an orthonormal matrix Q such that$\Vert (I-QQ^*)A\Vert\leq \epsilon$, with probability at least $1 − \min\{m, n\}10^{-r}.$
\item \textbf{randPI.m} - This is implementation of Algorithm 4.3 from paper \cite{4}.
		
		
		This function takes an $m \times n$ matrix A and integers l and q, it computes an $m \times l$ orthonormal matrix Q whose range approximates the range of A.
\item \textbf{randSI.m} -  This is implementation of Algorithm 4.4 from paper \cite{4}.
		
		
This function takes an $m \times n$ matrix A and integers l and q, it computes an $m \times l$ orthonormal matrix Q whose range approximates the range of A.
\item \textbf{FastRandRF.m} -  This is implementation of Algorithm 4.5 from paper \cite{4}.
		
		
		This function takes an $m \times n$ matrix A and an integer l, it computes an $m \times l$ orthonormal matrix Q whose range approximates the range of A. This function involved with sub sampled random Fourier transform (SRFT), see function \textbf{SRFT.m} in the Sub-Function session. 
\item \textbf{DirectEigvalueDecopo.m} - This is implementation of Algorithm 5.3 from paper \cite{4}.
		
		This function takes an Hermitian matrix A and a basis Q that can be generated by \textbf{FixedRank.m}, this computes an approximate eigenvalue decomposition $A \approx U\Lambda U^∗$, where U is orthonormal, and $\Lambda$ is a real diagonal matrix.
		
		To generate an input hermitian matrix, please visit and download the \href{https://www.mathworks.com/matlabcentral/fileexchange/25912-random-hermitian-matrix-generator}{\color{blue} \myul[blue] {random hermitian matrix generator function}}\cite{rherm}.
\item \textbf{EigvalueDecopoRow.m} - This is implementation of Algorithm 5.4 from paper \cite{4}.

        This function takes an Hermitian matrix A and a basis Q that can be generated by {FixedRank.m}, this computes an approximate eigenvalue decomposition  A≈UΛU∗ , where U is orthonormal, and  $\Lambda$  is a real diagonal matrix. This function is faster than \textbf{DirectEigvalueDecopo.m} but less accurate.
        
        To generate an input hermitian matrix, please visit and download the \href{https://www.mathworks.com/matlabcentral/fileexchange/25912-random-hermitian-matrix-generator}{\color{blue} \myul[blue] {random hermitian matrix generator function}}\cite{rherm}.
\item \textbf{EigvalueDecopoNystrom.m} - This is implementation of Algorithm 5.5 from paper \cite{4}.

        This function takes a positive semidefinite matrix A and a basis Q that can be generated by {FixedRank.m}, this computes an approximate eigenvalue decomposition  A≈UΛU∗ , where U is orthonormal, and  $\Lambda$  is nonnegative and diagonal.
\item \textbf{EigvalueDecopoOnePass.m} - This is implementation of Algorithm 5.6 from paper \cite{4}.

        This function takes an Hermitian  matrix A, a random test matrix $\Omega$, a sample matrix $Y =
A\Omega$, and an orthonormal matrix Q that can be generated by {FixedRank.m}, this computes an approximate eigenvalue decomposition  A≈UΛU∗. This algorithm requires only one pass for the input matrix A; Comparing to previous algorithms it takes less storage.  
\item \textbf{BasicRandSVD.m} - This in implementation of RSVD from paper \cite{5}.
		
      This function takes an $m \times n$ matrix A, a target rank k, and an over-sampling parameter p and computes matrices U, D, and V in an approximate rank-$(k + p)$ SVD of A (so that U and V are
orthonormal, D is diagonal, and $A\approx UDV^* $

    To generate an input hermitian matrix, please visit and download the \href{https://www.mathworks.com/matlabcentral/fileexchange/25912-random-hermitian-matrix-generator}{\color{blue} \myul[blue] {random hermitian matrix generator function}}\cite{rherm}.
\item \textbf{AERandSVD.m} - This is implementation of ALGORITHM: ACCURACY ENHANCED RANDOMIZED SVD from paper \cite{5}.
		
      This function takes an $m \times n$ matrix A, a target rank k, an over-sampling parameter p, and a exponent q. It computes matrices U, D, and V in an approximate rank-$(k + p)$ SVD of A (so that U and V are orthonormal, D is diagonal, and $A\approx UDV^* $. This algorithm is more accurate Compares to \textbf{BasicRandSVD.m}.
\item \textbf{AEORandSVD.m} - This is implementation of ALGORITHM: ACCURACY ENHANCED RANDOMIZED SVD (WITH ORTHONORMALIZATION) from paper \cite{5}.
		
      This function takes an $m \times n$ matrix A, a target rank k, an over-sampling parameter p, and a exponent q. It computes matrices U, D, and V in an approximate rank-$(k + p)$ SVD of A (so that U and V are orthonormal, D is diagonal, and $A\approx UDV^* $. This algorithm reduces the truncating error from the power iteration with large q in \textbf{AERandSVD.m}.
\item \textbf{SPRandEVDH.m} - This is implementation of ALGORITHM: SINGLE-PASS RANDOMIZED EVD FOR A HERMITIAN MATRIX from paper \cite{5}.
		
      This function takes an $n \times n$ Hermitian matrix A, a target rank k, and an over-sampling parameter p and computes matrices U and D in an approximate rank-k EVD of A (so that UU is an orthonormal matrix, D is a diagonal matrix, and $A\approx UDU^* $. This algorithm requires only one pass or the input matrix A; Comparing to previous algorithms it takes less storage.
      
      To generate an input hermitian matrix, please visit and download the \href{https://www.mathworks.com/matlabcentral/fileexchange/25912-random-hermitian-matrix-generator}{\color{blue} \myul[blue] {random hermitian matrix generator function}}\cite{rherm}.
\item \textbf{SPRandSVD.m }- This is implementation of ALGORITHM: SINGLE-PASS RANDOMIZED SVD FOR A GENERAL MATRIX from paper \cite{5}.
		
       This function takes an $m \times n$  matrix A, a target rank k, and an over-sampling parameter p and computes matrices U, D, and V in an approximate rank-k SVD of A (so that U and V are orthonormal, D is diagonal, and $A\approx UDV^* $. This function extends \textbf{SPRandEVDH.m} to general input matrix.
       
\item \textbf{randPowerMethod.m} -  This is implementation of Algorithm 4 from paper \cite{6}.		
    
			This function takes a Hermitian matrix A, a  number q for maximum number of iterations and a stopping tolerance $\epsilon$, it computes estimated $\xi$ for a maximum eigenvalue of A.
			
			To generate an input hermitian matrix, please visit and download the \href{https://www.mathworks.com/matlabcentral/fileexchange/25912-random-hermitian-matrix-generator}{\color{blue} \myul[blue] {random hermitian matrix generator function}}\cite{rherm}.
\item \textbf{randomizedLanczos.m} -  This is implementation of Algorithm 5 from paper \cite{6}.		
    
			This function takes a Hermitian matrix A, a  number q for maximum number of iterations and computes estimated $(\xi; y)$ for a maximum eigenpair of A.
            
            To generate an input hermitian matrix, please visit and download the \href{https://www.mathworks.com/matlabcentral/fileexchange/25912-random-hermitian-matrix-generator}{\color{blue} \myul[blue] {random hermitian matrix generator function}}\cite{rherm}.
\end{itemize}

\begin{table}[H]
\begin{tabular}{ | c | c | c |}
\hline
	  \thead{Function name} &   \thead{Reference}  & \thead{Algorithm \# or name\\ listed in the reference} \\ \hline
	EigvalueDecopoRow & \cite{2} & 5.4 \\ \hline
	EigvalueDecopoOnePass & \cite{2} & 5.6 \\ \hline
	EigvalueDecopoNystrom & \cite{2} & 5.5 \\ \hline
	DirectEigvalueDecopo & \cite{2} & 5.3 \\ \hline
	FastRandRF & \cite{2} & 4.5 \\ \hline
	Adaptive\_randRF & \cite{2} & 4.2 \\ \hline
	\makecell{AEORandSVD} & \cite{5} & \makecell{ACCURACY ENHANCED RANDOMIZED SVD \\(WITH ORTHONORMALIZATION) } \\ \hline
	AERandSVD & \cite{5} & \makecell{ACCURACY ENHANCED RANDOMIZED SVD} \\ \hline
	BasicRandSVD & \cite{5} & RSVD \\ \hline
	FixedRank & \cite{2} & \makecell{Proto-Algorithm: Solving the Fixed\\ -Rank Problem} \\ \hline
	randomizedLanczos & \cite{6} & 5 \\ \hline
	randPI & \cite{2} & 4.3 \\ \hline
	randPowerMethod & \cite{6} & 4 \\ \hline
	randRF & \cite{2} & 4.1 \\ \hline
	randSI & \cite{2} & 4.4 \\ \hline
	RandSVD & \cite{2} & Prototype for Randomized SVD \\ \hline
	SPRandEVDH & \cite{5} & \makecell{SINGLE-PASS RANDOMIZED EVD \\FOR A HERMITIAN MATRIX} \\ \hline
	SPRandSVD & \cite{5} & \makecell{SINGLE-PASS RANDOMIZED SVD \\FOR A GENERAL MATRIX} \\ \hline

\end{tabular}
\end{table}
\paragraph{Sub Functions:}
\begin{itemize}
\item \textbf{SRFT.m}  -  This function proceed the subsampled random Fourier transform and produce an $n\times l$ matrix $\Omega$ in the form \[\Omega = \sqrt{\frac{n}{l}}DFR,\]where 
\begin{itemize}
\begin{itemize}
\item D is an $n\times n$ diagonal matrix whose entries are independent random variables
uniformly distributed on the complex unit circle.
\item F is the $n\times n$ unitary discrete Fourier transform, whose entries take the values $f_{pq} = n^{-\frac{1}{2}}e^{-2\pi i(p-1)(q-1)/n}$ for $ p,q = 1,2,...,n$.
\item R is an $n\times l$ matrix that samples l coordinates from n uniformly at random.

\end{itemize}
\end{itemize}
\item \textbf{SPD.m} - This function randomly generates a symmetric positive defined matrix A with dimension $n\times n$.
\item \textbf{maxeig.m} - This function returns the maximum absolute eigenvalue of input matrix A.

\end{itemize}


\begin{table}[H]
\begin{tabular}{ | c | c | }
\hline
	\textbf{Function name} & \textbf{For} \\ \hline
	SRFT & subsampled random Fourier transform \\ \hline
	SPD & generates a random spd matrix \\ \hline
	maxeig & returns the max eigenvalue \\ \hline
\end{tabular}
\end{table}



\paragraph{Drivers:}
\begin{itemize}
\item \textbf{imagedriver.m} - Image reconstruction using different methods of randomized SVD. This driver produce all experiment results for experiment 4.1 \textit{Illustration of various Randomized SVD} listed in this report. 
\item \textbf{driver\_bound.m} - Analysis on error bounds of different methods of the 1st stage of randomized SVD. This driver produce all experiment results for experiment 4.2 \textit{Illustration of Error Bounds} listed in this report. 

To compile and run this test driver, please download the foloowing file from GitHub:
\begin{itemize}
\item \href{https://github.com/arvindks/randsvs/blob/master/testmatrices/controlledgap.m
}{\color{blue} \myul[blue] {controlledgap.m}} - This function takes input: $m,n$ as the size of desired matrix A, $r$ as the position of the gap and $gap$ for the size of the gap, and returns the testing matrix A that contains such a gap between its singular values at the defined position. Introduced in paper \cite{saibaba2019randomized}.

\end{itemize}
\item \textbf{driver\_sin.m} - Canonical angles for different methods of randomized SVD. This driver produce all experiment results for experiment 4.3 \textit{Illustration of Canonical angles} listed in this report. 

To compile and run this test driver, please download the following files from GitHub:
\begin{itemize}

\item \href{https://github.com/arvindks/randsvs/blob/master/testmatrices/controlledgap.m
}{\color{blue} \myul[blue] {controlledgap.m}} - This function takes input: $m,n$ as the size of desired matrix A, $r$ as the position of the gap and $gap$ for the size of the gap, and returns the testing matrix A that contains such a gap between its singular values at the defined position. Introduced in paper \cite{saibaba2019randomized}.
\item \href{https://github.com/arvindks/randsvs/blob/master/core/angle_bounds.m}{\color{blue} \myul[blue] {angle\_bounds.m}} - This function takes the right singular vectors V, the starting guess $\Omega$, the singular values s, a target rank k and the number of subspaces q, outputs both the bounds for sin\((\theta(U_1,U_h))\) and sin\((\theta(V_1,V_h))\). Introduced in paper \cite{saibaba2019randomized}.
\item \href{https://github.com/arvindks/randsvs/blob/master/core/subspace_angles.m}{\color{blue} \myul[blue] {subspace\_angles.m}}  - This function computes the canonical angles between two subspaces $U$ and $U_h$. Introduced in paper\cite{saibaba2019randomized}.
\end{itemize}
\end{itemize}
\begin{table}[H]
\begin{tabular}{ | c | c | c |}
\hline
	\textbf{Driver file name} & \textbf{For}  \\ \hline
	imagedriver & image experiment on chapter 3 \\ \hline
	driver\_bound & bound experiment on chapter 3  \\ \hline
	driver\_sin & canonical angles experiment on chapter 3 \\ \hline
\end{tabular}
\end{table}


\begin{table}[H]
\begin{tabular}{ | c | c | c |}
\hline
	\textbf{Driver file name} & \textbf{For}  \\ \hline
	imagedriver & image experiment on chapter 3 \\ \hline
	driver\_bound & bound experiment on chapter 3  \\ \hline
	driver\_sin & canonical angles experiment on chapter 3 \\ \hline
\end{tabular}
\end{table}

\paragraph{Data files:}
\begin{itemize}
\item \textbf{Sunflower.txt} - An $804\times 1092$ matrix converted from a photo of Kansas Sunflowers. 
\end{itemize}
\begin{table}[H]
\begin{tabular}{ | c | c | }
\hline
	\textbf{File name} & \textbf{For} \\ \hline
	Sunflower & matrix converted from a photo of Kansas Sunflowers \\ \hline
\end{tabular}
\end{table}

