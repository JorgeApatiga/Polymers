library(sqldf) 
Tiempo_percolacion<-function(x1,x2,xx,y,z,w,xy){ 
  #Codigo para generar la red de Grafeno 
  
  #Primero la dimencion de la red es nxm 
  n<-z 
  m<-w 
  theta1<-x1 
  theta2<-x2 
  veci<-function(u){   
    if(((u[1]-1)<1) & ((u[2]-1)<1)){  #esqina inferior izquierda 
      r1<-c(u[1]+1,u[2])   
      r2<-c(u[1],u[2]+1) 
      r3<-c(u[1]+1,u[2]+1) 
      vecinos<-rbind(r1,r2,r3) 
      return(vecinos) 
    }else if(((u[1]+1)>n) & ((u[2]+1)>m)){  #esquina superior derecha 
      r1<-c(u[1]-1,u[2])   
      r2<-c(u[1],u[2]-1) 
      r3<-c(u[1]-1,u[2]-1) 
      vecinos<-rbind(r1,r2,r3) 
      return(vecinos) 
    }else if(((u[1]-1)<1) & ((u[2]+1)>m)){  #esquina superior izquierda 
      r1<-c(u[1]+1,u[2])   
      r2<-c(u[1],u[2]-1) 
      r3<-c(u[1]+1,u[2]-1) 
      vecinos<-rbind(r1,r2,r3) 
      return(vecinos) 
    }else if(((u[1]+1)>n) & ((u[2]-1)<1)){  #esquina inferior derecha 
      r1<-c(u[1]-1,u[2])   
      r2<-c(u[1],u[2]+1) 
      r3<-c(u[1]-1,u[2]+1) 
      vecinos<-rbind(r1,r2,r3) 
      return(vecinos) 
    }else if((u[1]-1)<1){  #primera columna 
      r1<-c(u[1]+1,u[2]) 
      r2<-c(u[1],u[2]+1) 
      r3<-c(u[1]+1,u[2]+1) 
      r4<-c(u[1],u[2]-1) 
      r5<-c(u[1]+1,u[2]-1) 
      vecinos<-rbind(r1,r2,r3,r4,r5) 
      return(vecinos) 
    }else if((u[2]-1)<1){ #primera fila 
      r1<-c(u[1]-1,u[2]) 
      r2<-c(u[1]+1,u[2]) 
      r3<-c(u[1]+1,u[2]+1) 
      r4<-c(u[1],u[2]+1) 
      r5<-c(u[1]-1,u[2]+1 ) 
      vecinos<-rbind(r1,r2,r3,r4,r5) 
      return(vecinos) 
    }else if((u[1]+1)>n){  #ultima columna 
      r1<-c(u[1]-1,u[2]) 
      r2<-c(u[1]-1,u[2]-1) 
      r3<-c(u[1]-1,u[2]+1) 
      r4<-c(u[1],u[2]+1) 
      r5<-c(u[1],u[2]-1 ) 
      vecinos<-rbind(r1,r2,r3,r4,r5) 
      return(vecinos) 
    }else if((u[2]+1)>m){  #ultia fila 
      r1<-c(u[1]+1,u[2]) 
      r2<-c(u[1]-1,u[2]) 
      r3<-c(u[1]-1,u[2]-1) 
      r4<-c(u[1],u[2]-1) 
      r5<-c(u[1]+1,u[2]-1) 
      vecinos<-rbind(r1,r2,r3,r4,r5) 
      return(vecinos) 
    }else{                #puntos centrales 
      r1<-c(u[1]-1,u[2]) 
      r2<-c(u[1]-1,u[2]+1) 
      r3<-c(u[1],u[2]+1) 
      r4<-c(u[1]+1,u[2]+1) 
      r5<-c(u[1]+1,u[2]) 
      r6<-c(u[1]+1,u[2]-1) 
      r7<-c(u[1],u[2]-1) 
      r8<-c(u[1]-1,u[2]-1) 
      vecinos<-rbind(r1,r2,r3,r4,r5,r6,r7,r8) 
      return(vecinos)  
    } 
  } 
  
  Dat<-rep(0,((n+xx)*m)) 
  dat<-matrix(Dat, nrow=(n+xx),ncol=m,byrow=FALSE) 
  for(i in 1:(n+xx)){ 
    for(j in 1:m){ 
      aux1<-veci(c(i,j)) 
      aux2<-0 
      for(ii in 1:length(aux1[,1])){ 
        if(dat[aux1[ii,1],aux1[ii,2]]!=0){ 
          aux2<-aux2+1 
        } 
      } 
      if(aux2!=0){dat[i,j]<-rbinom(1,1,theta1)}else{dat[i,j]<-rbinom(1,1,theta2)} 
    } 
  } 
  
  dat<-dat[-(1:xx),] 
  
  #Primero un vector auxiliar donde guardamos las entradas de los clusters 
  cl<-c(0,0) 
  k<-0  #numerador de numero de clusters 
  
  #Comensam os el ciclo que recorre el polimero y nos regresa los clusters. 
  for(i in 1:n){ 
    for(j in 1:m){ 
      if(dat[i,j]==1){   #queremos que esa entrada no sea 0 y no este en un cluster anterior 
        pre_clus<-matrix(c(i,j),ncol=2,byrow = TRUE)              # vector con las coordenadas que estamos operando 
        pre_clus_save<-matrix(c(0,0),ncol=2,byrow = TRUE) 
        vecinos_1<-matrix(c(i,j),ncol=2,byrow = TRUE)             #matriz auxiliar que guarda 1 
        aux1<-data.frame(pre_clus)     #los convertimos en data.frame para sacar la diferencia y regresamos matries 
        aux2<-data.frame(pre_clus_save)    # HAY QUE INSTALAR EL PAQUETE SQLDF###### 
        matrix.diff<-data.matrix(sqldf("select * from [aux1] EXCEPT SELECT * FROM [aux2]")) 
        while(length(matrix.diff[,1])>0){        
          pre_clus_save<-pre_clus 
          for(l in 1:length(matrix.diff[,1])){ 
            aux3<-veci(c(matrix.diff[l,1],matrix.diff[l,2])) 
            for(p in 1:length(aux3[,1])){      
              aux4<-c(aux3[p,1],aux3[p,2]) 
              if(dat[aux4[1],aux4[2]]==1){ 
                vecinos_1<-rbind(vecinos_1,aux4) 
              } 
            } #nos regresa un vector pre_clus con los vecinos que tienen valor 1 
            aux5<-data.frame(vecinos_1)       #convertimos en data frame para diferenciar 
            aux6<-data.frame(pre_clus) 
            nuevos_puntos<-data.matrix(sqldf("select * from [aux5] EXCEPT SELECT * FROM [aux6]")) 
            if(length(nuevos_puntos[,1])>0){pre_clus<-rbind(pre_clus,nuevos_puntos)} #actualiza nuestro vector princial 
          } 
          aux1<-data.frame(pre_clus)     #si en el if anterior no habia nuevos vecinos con 1 se queda igual la matriz 
          aux2<-data.frame(pre_clus_save)    #como ya actualizamos abajo del while volvemos a comparar 
          matrix.diff<-data.matrix(sqldf("select * from [aux1] EXCEPT SELECT * FROM [aux2]"))  
        } #ciclo mientras hay vecinos con unos(1) 
        if(length(pre_clus[,1])>=y){            #si el numero de elementos en pr_clus es mayor a 5 es un cluster con numero k 
          k<-k+1 
          for(v in 1:length(pre_clus[,1])){ 
            dat[pre_clus[v,1],pre_clus[v,2]]<-(k+1) 
          } #cambia el valor de los elentos del cluster de 1 a (k+1) 
        }else{ 
          for(v in 1:length(pre_clus[,1])){dat[pre_clus[v,1],pre_clus[v,2]]<-0}  
        } #cambia el valor de los elementos de pre_clus a 0   
      }   #si el valor es 1 y el numero de elementos del cluster al que pertence es mayor a 5 los cambia  
      #a (k+1) de lo contrario a 0 
      #print(c(i,j))                        #########Punto de prueba 
    }  
  }   
  
  #saturacion 
  sat<-(sum(dat!=0)/(n*m))*100 
  sat 
  
  
  #codigo auxiliar para los colores 
  colores<-function(x){ 
    if(x==0){return("gray")} 
    else{return("black")} 
  } 
  
  #Extraemos las coordenadas de los puntos para graficarlos 
  plot(1,1,main = paste("Depositos de Grafeno",perc),xlim=c(1,n), ylim =c(1,m)) 
  for(i in 1:n){#Esto lo hago para que grafique cuadros 
    for(j in 1: m){ 
      rect(i-1,j-1,i,j, col=colores(dat[i,j]),border = NA) # transparent  
    } 
  } 
  
  
  
  ####Hacer cumulos de 3 particulas y umbral de 4 espacios. 
  Mf<-list()               ### Mf es una lista que contiene los clusters. 
  
  
  #primero haremos matriz de cada cluster y una lista de matrices para comparar. 
  for(p in 1:k){ 
    M<-matrix(c(0,0),nrow=1) 
    nam<-paste("M", p, sep = "") 
    for(j in 1:m){ 
      for(i in 1:n){ 
        if(dat[i,j]==(p+1)){ 
          aux7<-c(i,j) 
          M<-rbind(M,aux7) 
        }   
      } 
    } 
    M<-M[-1,] 
    Mf[[p]]<-M 
    assign(nam, M) 
  } 
  
  #Ahora  calculamos la distancia entre las matrices. 
  eucd <- function(p,q){ 
    sqrt(sum((p - q)^2)) 
  } 
  
  Mdivsw<-NULL 
  Mdist<-matrix(c(0,0,0), ncol=3, byrow = TRUE) 
  MDIST <- matrix(data=NA,nrow=(k),ncol=(k)) 
  
  for( i in 1:(k-1)){ 
    matrix1ength1<-length(Mf[[i]][,1]) 
    for(w in (i+1):k){ 
      Mdivsw<-c(1000) 
      matrixlength2<-length(Mf[[w]][,1]) 
      for (j in 1:matrix1ength1) { 
        for(q in 1:matrixlength2){ 
          dist<-eucd((Mf[[i]][j,]),(Mf[[(w)]][q,])) 
          Mdivsw<-c(Mdivsw,dist) 
        } 
      } 
      mindis<-min(Mdivsw) 
      Mdist<-rbind(Mdist,c(i,w,mindis)) 
      MDIST[(i),(w)]<-mindis 
      MDIST[(w),(i)]<-mindis 
    } 
  } 
  
  
  # Ahora determinaremos los clusters que estan a una distancia no mas de 3 del inicio y no mas de 3 del final 
  #recordemos que el número de clusters es K 
  # El vector que creamos en este paso nos dice si el cluster p esta al inicio del polimero 
  clus_inicial<-vector() 
  clus_inicial_num<-vector() 
  
  for(p in 1:k){ 
    if(sum(Mf[[p]][,1]<xy)>0){ 
      clus_inicial<-c(clus_inicial,1) 
      clus_inicial_num<-c(clus_inicial_num,p) 
    }else{ 
      clus_inicial<-c(clus_inicial,0) 
    } 
  } 
  clus_inicial 
  if(sum(clus_inicial)==0){ 
    return(c(sat,0)) 
    break 
  }   #si no hay cluster inicial sales 
  
  # Ahora creamos el vector de probabilidad inicial de manera uniforme 
  
  Pi_inicial<-clus_inicial/(sum(clus_inicial)) 
  
  
  # Ahora  el vector con clusters al final 
  
  tau_pre_absorvente<-vector() 
  tau_pre_absorvente_num<-vector() 
  
  for(p in 1:k){ 
    if(sum(Mf[[p]][,1]>(n-xy))>0){ 
      tau_pre_absorvente<-c(tau_pre_absorvente,p) 
      tau_pre_absorvente_num<-c(tau_pre_absorvente_num,1) 
    }else{ 
      tau_pre_absorvente_num<-c(tau_pre_absorvente_num,0) 
    } 
  } 
  tau_pre_absorvente 
  if(sum(tau_pre_absorvente_num)==0){ 
    return(c(sat,0)) 
    break 
  }    #si no hay estado absorbente no conduce y sale. 
  
  
  
  #creamos una matris para poder generar nuestra matriz de intencidad 
  aux8<-rep(0,(k*k)) 
  Mat_intensidad_pre_absorbente<-matrix(aux8,nrow = k,ncol = k) 
  
  for(j in 1:k){ 
    for(i in 1:k){ 
      if(i < j){ 
        Mat_intensidad_pre_absorbente[i,j]<-(1/(MDIST[i,j]^(12)))*(MDIST[i,j]>(sqrt(n^2+m^2))*0.3)+(1/(MDIST[i,j]^6))*(MDIST[i,j]<=(sqrt(n^2+m^2))*0.3) 
      }else if(i > j){ 
        Mat_intensidad_pre_absorbente[i,j]<-((1/(MDIST[i,j]^12))*(MDIST[i,j]>(sqrt(n^2+m^2))*0.3)+(1/(MDIST[i,j]^6)))*(MDIST[i,j]<=(sqrt(n^2+m^2))*0.3) 
      }else{ 
        Mat_intensidad_pre_absorbente[i,j]<-0 
      } 
    } 
  } 
  
  #asignamos a la diagonal como menos la suma del resto de las columnas de cada renglon 
  Mat_intensidad_pre_absorbente_diagonal<-vector() 
  for(i in 1:k){ 
    Mat_intensidad_pre_absorbente_diagonal<-c(Mat_intensidad_pre_absorbente_diagonal,-sum(Mat_intensidad_pre_absorbente[i,])) 
  } 
  for(i in 1:k){ 
    Mat_intensidad_pre_absorbente[i,i]<-Mat_intensidad_pre_absorbente_diagonal[i] 
  } 
  
  
  
  #En  Mat_subintensidad tenemos la matriz sin los renglones del estado absorbente. 
  Mat_subintensidad<-Mat_intensidad_pre_absorbente[-tau_pre_absorvente,]  
  #Sumamos las probabilidades de los estados absorbentes y los mandamos a tau 
  tau_absorvente<-vector() 
  for(i in 1:length(Mat_subintensidad[,1])){ 
    tau_absorvente<-c(tau_absorvente,sum(Mat_subintensidad[i,tau_pre_absorvente])) 
  } 
  #eliminamos las columnas del estado absorbente de la matris Mat_intensidad_absorbente 
  Mat_subintensidad<-Mat_subintensidad[,-tau_pre_absorvente] 
  
  #Finalmente creamos la matriz de intensidad. 
  tau<-matrix(tau_absorvente,ncol=1) 
  Mat_intensidad<-cbind(Mat_subintensidad,tau) 
  ren_final<-rep(0,length(Mat_intensidad[1,])) 
  Mat_intensidad<-rbind(Mat_intensidad,ren_final) 
  if(rcond(-Mat_subintensidad) < .Machine$double.eps){ 
    return(c(sat,0)) 
    break 
  }  
  
  U_Mat<-solve(-Mat_subintensidad)  
  tiempo_final<-vector() 
  for(i in 1:length(clus_inicial_num)){ 
    tiempo_final<-c(tiempo_final,sum(U_Mat[clus_inicial_num[i],])) 
  } 
  Tiempo_final<-sum(tiempo_final/sum(clus_inicial)) 
  return(c(sat,Tiempo_final)) 
} 

prom<-function(x){ 
  j<-1 
  jl<-0 
  Tabla_promedios<-matrix(c(0,0,0,0),ncol=4) 
  jj<-(round(x[length(x[,1]),1])-round(x[2,1]))*5+1 
  for(i in 1:jj){ 
    aux1<-vector() 
    aux2<-vector() 
    time<-(round(x[2,1]))+(i*0.2) 
    while ((time>x[j,1]) & x[j,1]<=x[length(x[,1]),1]){ 
      if(x[j,2]!=0){ 
        aux1<-c(aux1,x[j,2]) 
      } 
      if((j+1)<=length(x[,1])){ 
        j<-j+1 
      }else{ 
        break 
      } 
    } 
    aux2<-c((time-0.2),time,1/mean(aux1),j-jl) 
    Tabla_promedios<-rbind(Tabla_promedios,aux2) 
    jl<-j 
  } 
  return(Tabla_promedios[-1,]) 
} 
Tabla_final<-matrix(c(0,0),ncol = 2,byrow = TRUE) 
for(perc in 1:10000){ 
  Tabla_final<-rbind(Tabla_final,Tiempo_percolacion(0.38,0.34,20,3,100,20,6)) 
  print(perc) 
  if(perc%%1000==0){ 
    write.table(Tabla_final,file=" C:/Users/jl_ap/OneDrive/Documentos/Articulos/Roxana/Grafeno/Programa/Resultados /Grapheno_1.txt",row.names=F, sep="\t") 
    Tabla_prom<-prom(Tabla_final) 
    write.table(Tabla_prom,file=" C:/Users/jl_ap/OneDrive/Documentos/Articulos/Roxana/Grafeno/Programa/Resultados /Grapheno_1_prom.txt",row.names=F, sep="\t") 
  } 
} 
Tabla_final<-Tabla_final[order(Tabla_final[,1]),] 
Tabla_final<-Tabla_final[-1,] 
write.table(Tabla_final,file=" C:/Users/jl_ap/OneDrive/Documentos/Articulos/Roxana/Grafeno/Programa/Resultados /Grapheno_1.txt",row.names=F, sep="\t") 



Tabla_prom<-prom(Tabla_final) 
Tabla_prom 
write.table(Tabla_prom,file=" C:/Users/jl_ap/OneDrive/Documentos/Articulos/Roxana/Grafeno/Programa/Resultados /Grapheno_1_prom.txt",row.names=F, sep="\t") 
