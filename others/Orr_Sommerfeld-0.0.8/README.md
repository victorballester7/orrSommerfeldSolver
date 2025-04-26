#Cmd Packaging libs in python

# to upgrade tools 
> python3 -m pip install --user --upgrade setuptools wheel
> python3 -m pip install --user --upgrade twine


#to upload library:
>cd $folder_lib_name

>python3 setup.py sdist bdist_wheel	#Genera los archivos para subir a www.pypi.org
>twine upload dist/* 			#Sube los archivos para distribución, previo configurar una cuenta en www.pypi.org 


#to install de library
pip install $lib_name==$last_version			#se descarga en este caso se creo la demo para el proy integrador: calc_demo_pi==0.05


#using lib: (eg: ipython)
In [1]: from calc_demo_pi.Chebychev import func    #se importa la bibl 
In [2]: func.run()				   #corre una funcion print.


Out[2]: 'Successful!'

#-----------------------------------------

In [2]: func.cheb(5)                                                            
Out[2]: 
(array([[  8.5       , -10.47213595,   2.89442719,  -1.52786405,
           1.10557281,  -0.5       ],
        [  2.61803399,  -1.17082039,  -2.        ,   0.89442719,
          -0.61803399,   0.2763932 ],
        [ -0.7236068 ,   2.        ,  -0.17082039,  -1.61803399,
           0.89442719,  -0.38196601],
        [  0.38196601,  -0.89442719,   1.61803399,   0.17082039,
          -2.        ,   0.7236068 ],
        [ -0.2763932 ,   0.61803399,  -0.89442719,   2.        ,
           1.17082039,  -2.61803399],
        [  0.5       ,  -1.10557281,   1.52786405,  -2.89442719,
          10.47213595,  -8.5       ]]),
 array([ 1.        ,  0.80901699,  0.30901699, -0.30901699, -0.80901699,
        -1.        ]))
