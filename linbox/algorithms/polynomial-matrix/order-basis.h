/* linbox/algorithms/sigma-basis.h
 * Copyright (C) 2005,2013 Pascal Giorgi, Romain Lebreton
 *
 * Written by Pascal Giorgi pascal.giorgi@lirmm.fr
 *            Romain Lebreton lebreton@lirmm.fr
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */





#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/polynomial-matrix.h"

#ifdef TRACK_MEMORY_MATPOL
#define MEMINFO2 STR_MEMINFO<<MEMINFO
#else
#define MEMINFO2 ""
#endif
#include "linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h"
#include <vector>
#include <algorithm>
#include <fstream>
#include <chrono>
#include "fflas-ffpack/fflas-ffpack.h"
#define MBASIS_THRESHOLD_LOG 5
#define MBASIS_THRESHOLD (1<<MBASIS_THRESHOLD_LOG)



namespace LinBox {


#ifdef __MINPOLY_SETTING
#define __CHECK_PMBASIS
#define __PROBA_CHECK
#define __DUMP_ORDERBASIS
#define __CHECK_PMBASIS_THRESHOLD 1023
#define  PROFILE_PMBASIS
#endif

        
#ifdef __CHECK_ORDERBASIS
#define __CHECK_MBASIS
#define __CHECK_PMBASIS
#endif

#ifndef __CHECK_PMBASIS_THRESHOLD 
#define __CHECK_PMBASIS_THRESHOLD MBASIS_THRESHOLD
#endif        

#if defined (__CHECK_MBASIS) or defined (__CHECK_PMBASIS)
#include <string>
        template<typename Field, typename Mat>
        std::string check_orderbasis(const Field& F, const Mat& sigma,  const Mat& serie, size_t ord, int val,std::vector<size_t> &shift){
                PolynomialMatrixMulDomain<Field> PMD(F);
                MatrixDomain<Field> MD(F);
#ifdef __PROBA_CHECK
                Mat T(F,sigma.rowdim(),1,sigma.size()+serie.size()-1);
                Mat U(F,serie.coldim(),1,1), serieU(F,serie.rowdim(),1,sigma.size()+serie.size()-1);
                typename Field::RandIter Gen(F);
                for (size_t i=0;i<serie.coldim();i++)
                        Gen.random(U.ref(i,0,0));
                PMD.mul(serieU,serie,U);                
                PMD.mul(T,sigma,serieU);
#else
                Mat T(F,sigma.rowdim(),serie.coldim(),sigma.size()+serie.size()-1);
                PMD.mul(T,sigma,serie);
#endif

                size_t i=0;
                std::string msg(".....");
                bool nul_sigma=true;
                while(i<ord && MD.isZero(T[i])){
                        if (i<sigma.size() && !MD.isZero(sigma[i])) nul_sigma=false;		
                        i++;
                }
                if (i<ord){
                        std::cout<<"error at degree="<<i<<std::endl;
                        T[i].write(std::cout, Tag::FileFormat::Plain);
                        std::cout<<"***"<<std::endl;
#ifdef __DEBUG_ORDERBASIS
                     
                        std::cout<<serie<<std::endl;
                        std::cout<<sigma<<std::endl;
#endif
#ifndef __DUMP_ORDERBASIS
                        std::terminate();
#endif
                }
	
	
                if (i==ord && !nul_sigma){
                        msg+="done";
#ifdef __DUMP_ORDERBASIS
                        std::string file_str("orderbasis-");
                        file_str+= std::to_string(ord)+"-"+std::to_string(val)+".dump";
                        std::ofstream file(file_str);
                        file<<"# order = "<<ord<<std::endl;
                        file<<"# shifted degree =";
                        std::copy(shift.begin(), shift.end(), std::ostream_iterator<size_t>(file, " "));
                        file<<"\n # orderbasis: \n";                        
                        sigma.dump(file);
                        file.close();
                        msg+="   ---> dumping it to "+file_str;                        
#endif
                }
                else {
                        msg+="error";
#ifdef __DUMP_ORDERBASIS
                        std::string file_str("orderbasis-bad-");
                        file_str+= std::to_string(ord)+"-"+std::to_string(val)+".dump";
                        std::ofstream file(file_str);
                        file<<"# order = "<<ord<<std::endl;
                        file<<"# shift =";
                        std::copy(shift.begin(), shift.end(), std::ostream_iterator<size_t>(file, " "));
                        file<<"\n # orderbasis: \n";                        
                        sigma.dump(file);
                        file.close();
                        msg+="   ---> dumping it to "+file_str;
                        std::cerr<<msg<<std::endl;
                        std::cerr<<"Aborting order basis computation ...\n";
                        std::terminate();
#endif
                }
                return msg;
        }
#endif

        
        template< size_t K>
        struct EarlyTerm {
                size_t _count;
                size_t _val;

                EarlyTerm():  _count(0),_val(0){}

                void update(size_t r, const std::vector<size_t>& u){
                        std::vector<size_t> v(u);
                        sort(v.begin(),v.end());
                        size_t x=0;
                        for (size_t i=0;i<r;i++)
                                x+=v[i];
                        if (x==_val)
                                _count++;
                        else{
                                _val=x;
                                _count=0;
                        }
                }

                bool terminated() const {return _count>=K;}

                void reset() {_count=0;_val=0;}
        };

		template<class Field, class ET=EarlyTerm<(size_t) -1> >
        class OrderBasis {
        public:
                typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,Field> MatrixP;
                typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field> PMatrix;
        private:
                const Field*                     _field;
                PolynomialMatrixMulDomain<Field>   _PMD;
                BlasMatrixDomain<Field>            _BMD;
                ET                           _EarlyStop;
        public:
#if  defined(PROFILE_PMBASIS) or defined(__CHECK_MBASIS) or defined(__CHECK_PMBASIS)
                size_t _idx=0;
                size_t _target=0;
                double  _eta=0.;
                std::chrono::time_point<std::chrono::system_clock> _start, _end;
                bool _started=false;
#endif
                OrderBasis(const Field& f) : _field(&f), _PMD(f), _BMD(f) {                 
                }

                inline const Field& field() const {return *_field;}

                // serie must have exactly order elements (i.e. its degree = order-1)
                // sigma can have at most order+1 elements (i.e. its degree = order)
                template<typename PMatrix1, typename PMatrix2>
                size_t PM_Basis2(PMatrix1                 &sigma,
                                 const PMatrix2           &serie,
                                 size_t                    order,
                                 std::vector<size_t>       &shift){
                        std::cout<<"COPYING INTIAL SERIE"<<std::endl;
                        PMatrix2 serie2(serie.field(),serie.rowdim(),serie.coldim(),serie.size());
                        serie2.copy(serie);
                        return PM_Basis(sigma,serie2,order,shift);
                }

                
                // serie must have exactly order elements (i.e. its degree = order-1)
                // sigma can have at most order+1 elements (i.e. its degree = order)
                template<typename PMatrix1, typename PMatrix2>
                size_t PM_Basis(PMatrix1                 &sigma,
                                const PMatrix2           &serie,
                                size_t                    order,
                                std::vector<size_t>       &shift)
                {

#ifdef PROFILE_PMBASIS
                        //std::cout<<"Start PM-Basis : "<<order<<" ("<<_idx<<"/"<<_target<<")] : "<<std::endl;//MEMINFO2<<std::endl;
                        if (_target==0) _target=order;
                        if (!_started) {_started=true; _start = std::chrono::system_clock::now();}
                        std::chrono::time_point<std::chrono::system_clock> _chrono_start=std::chrono::system_clock::now();
#endif
                        
                        if (order <= MBASIS_THRESHOLD) {
#if defined (PROFILE_PMBASIS) or defined(__CHECK_PMBASIS)
                                _idx+=order;
#endif
                                return M_Basis(sigma, serie, order, shift);                            
                        }
                        else {
#ifdef PROFILE_PMBASIS
                                Timer chrono;
                                chrono.start();
#endif
                                size_t ord1,ord2,d1,d2;
                                ord1 = order>>1;
                                ord2 = order-ord1; // ord1+ord2=order
                                //ord2 = order>>1;
                                //ord1 = order-ord2; // ord1+ord2=order
                                size_t m,n,k;
                                m=sigma.rowdim();
                                n=sigma.coldim();
                                k=serie.coldim();
                                integer p;

                                // first recursive call
                                PMatrix1 sigma1(field(),m,n,ord1+1);
                                
#ifdef MEM_PMBASIS
                                std::cerr<<"[PM-Basis ("<<order<<") "<<_idx<<"/"<<_target<<"] [Sigma1] -> "<<MB(sigma1.realmeminfo())<<"Mo"<<MEMINFO2<<std::endl;
#endif
                                PMatrix2 *serie1 = new PMatrix2(field(),n,k,ord1);
#ifdef MEM_PMBASIS
                                std::cerr<<"[PM-Basis ("<<order<<") "<<_idx<<"/"<<_target<<"] [Serie1] -> "<<MB(serie1->realmeminfo())<<"Mo"<<MEMINFO2<<std::endl;
#endif
                                serie1->copy(serie,0,ord1-1);
                                d1 = PM_Basis(sigma1, *serie1, ord1, shift);
                                //DEL_MEM(serie1->realmeminfo())
                                delete serie1;                                
                                if (_EarlyStop.terminated()){
                                        sigma=sigma1;
                                        return d1;
                                }

                                // compute the serie update
                                // TODO: for Block Wiedemann, this step can use only the first column of sigma
                                PMatrix2 *serie2=new PMatrix2(field(),n,k,ord2);//serie2 size=ord1+1 -> midproduct)
                                //ADD_MEM(serie2->realmeminfo());
#ifdef MEM_PMBASIS
                                std::cerr<<"[PM-Basis ("<<order<<") "<<_idx<<"/"<<_target<<"] [Serie2] -> "<<MB(serie2->realmeminfo())<<"Mo"<<MEMINFO2<<std::endl;
#endif              
                                _PMD.midproductgen(*serie2, sigma1, serie, true, ord1+1,ord1+ord2);
                                
#ifdef PROFILE_PMBASIS
                                //chrono.stop();
                                //std::cout<<"      -> serie update "<<sigma1.size()<<"x"<<order<<" --> "<<chrono.usertime()<<std::endl;//MEMINFO2<<std::endl;
                                //chrono.clear();chrono.start();
#endif
                                // second recursive call
                                
                                PMatrix1 sigma2(field(),m,n,ord2+1);
#ifdef MEM_PMBASIS
                                std::cerr<<"[PM-Basis("<<order<<") "<<_idx<<"/"<<_target<<"] [Sigma2] -> "<<MB(sigma2.realmeminfo())<<"Mo"<<MEMINFO2<<std::endl;
#endif
                                d2 = PM_Basis(sigma2, *serie2, ord2, shift);
                                delete serie2;                                 

                                // compute the result
                                _PMD.mul(sigma, sigma2, sigma1);
                                sigma.resize(d1+d2+1);                                
#ifdef PROFILE_PMBASIS
                                //chrono.stop();
                                //std::cout<<"      -> basis product "<<sigma1.size()<<"x"<<sigma2.size()<<" = "<<d1+d2+1<<" -->"<<chrono.usertime()<<MEMINFO2<<std::endl;
#endif

#ifdef __CHECK_PMBASIS
                                std::cerr<<"PMBASIS: order "<<order<<check_orderbasis(field(),sigma,serie,order,((_idx/order)+1)&1,shift)<<std::endl;
#endif
#ifdef PROFILE_PMBASIS
                                chrono.stop();
                                _end = std::chrono::system_clock::now();                                
                                std::chrono::duration<double> elapsed_beginning = _end-_start;
                                std::chrono::duration<double> elapsed_comp      = _end-_chrono_start;

                                double magicnumber=double(_target)/double(order)*log(double(_target)/double(order))/log(2.);
                                double tcomp = elapsed_comp.count();
                                double telap = elapsed_beginning.count();
                                
                                _eta=(_eta!=0.0?std::min(_eta,tcomp*magicnumber):tcomp*magicnumber);
                                std::cerr<<"[PM-Basis : "<<order<<" ("<<_idx<<"/"<<_target<<")] : "<<chrono.usertime()
                                         << " (ETA: "<< telap<<"s / "<<_eta<<"s)"<<MEMINFO2<<std::endl;
                                chrono.clear();chrono.start();
#endif


                                return d1+d2;
                        }
                }

                // serie must have exactly order elements (i.e. its degree = order-1)
                size_t M_Basis(MatrixP              &sigma,
                               const MatrixP        &serie,
                               size_t                 order,
                               std::vector<size_t>   &shift)
                {

                        PMatrix sigma1(field(),sigma.rowdim(),sigma.coldim(),order+1);
                        PMatrix serie1(field(),serie.rowdim(),serie.coldim(),order);
                        serie1.copy(serie,0,order-1);
                        size_t d= M_Basis(sigma1,serie1,order,shift);
                        sigma.resize(d+1);
                        sigma.copy(sigma1,0,d);
                        return d;

                }
                
                // serie must have exactly order elements (i.e. its degree = order-1)
                template<typename PMatrix1, typename PMatrix2>
                size_t M_Basis(PMatrix1              &sigma,
                               const PMatrix2        &serie,
                               size_t                 order,
                               std::vector<size_t>   &shift)
                {
#ifdef __DEBUG_MBASIS
                        std::cout<<"------------- mba : "<<order<<std::endl;
                        std::cout<<serie<<std::endl;
#endif
                        size_t m=serie.rowdim();
                        size_t n=serie.coldim();
                        size_t rank=0;
                        BlasMatrix<Field> delta(field(),m,n);
                        size_t max_degree=0;        // a bound on the row degree of sigma
                        std::vector<size_t> degree(m,0); // a bound on each row degree of sigma
                        //auto degree=shift;
                        //size_t max_degree=*std::max_element(degree.begin(),degree.end());

                        // set sigma to identity
                        for(size_t i=0;i<m*m;i++)
                                sigma.ref(i,0)=0;
                        for(size_t i=0;i<m;i++)
                                sigma.ref(i,i,0)=1;

                        BlasPermutation<size_t> Qt(m), P(n);
                        TransposedBlasMatrix<BlasPermutation<size_t> > Q(Qt);
                        typedef BlasSubmatrix<BlasMatrix<Field> > View;
                        size_t k=0;
                        for (; k<order && !_EarlyStop.terminated(); k++){
#ifdef __DEBUG_MBASIS
                                std::cout<<std::endl<<"****************** "<<k<<std::endl;
                                std::cout<<"shift=";std::copy(shift.begin(),shift.end(),std::ostream_iterator<size_t>(std::cout,","));std::cout<<std::endl;
                                std::cout<<"degree=";std::copy(degree.begin(),degree.end(),std::ostream_iterator<size_t>(std::cout,","));std::cout<<std::endl;                                
                                
#endif
                                // sort the shift in ascending order (-> minimize the shifted row-degree of sigma)
                                // -> store the permutation in Bperm
                                // -> permute the row degree at the same time
                                std::vector<size_t> perm(m);
                                for (size_t i=0;i<m;i++) perm[i]=i;
                                for (size_t i=0;i<m;++i) {
                                        size_t idx_min=i;
                                        for (size_t j=i+1;j<m;++j)
                                                if (shift[j]< shift[idx_min])
                                                        idx_min=j;
                                        std::swap( shift[i],  shift[idx_min]);
                                        std::swap(degree[i], degree[idx_min]);
                                        perm[i]=idx_min;
                                }
                                BlasPermutation<size_t> Bperm(perm);
#ifdef __DEBUG_MBASIS
                                std::cout<<"Bp=";
                                Bperm.write(std::cout,false);
                                std::cout<<std::endl;                                
                                std::cout<<"Bp.shift=";std::copy(shift.begin(),shift.end(),std::ostream_iterator<size_t>(std::cout,","));std::cout<<std::endl;
                                std::cout<<"Bp.degree=";std::copy(degree.begin(),degree.end(),std::ostream_iterator<size_t>(std::cout,","));std::cout<<std::endl;                                
#endif
                                // check if Qt is identity
                                //size_t lll=0;
                                //while(lll<m && Qt[lll]==lll) lll++;

                                // permute the row of the current sigma basis
                                // and compute the new discrepancy
                                //if (true){
                                //if (k==0 || lll!=m){
                                if (k==0 ){ //|| !Qt.isIdentity()){
                                        _BMD.mulin_right(Bperm, sigma[0]);
                                        _BMD.mul(delta,sigma[0],serie[k]);
                                        for(size_t i=1;i<=std::min(k,max_degree);i++){
                                                _BMD.mulin_right(Bperm, sigma[i]);
                                                _BMD.axpyin(delta,sigma[i],serie[k-i]);
                                        }
                                }
                                else{
                                        // if Qt is identity then the first rank rows of delta remain unchanged
                                        // the first rank rows of Qt.delta remain unchanged
                                        if (!Qt.isIdentity())
                                                _BMD.mulin_right(Qt, delta);

                                        View delta1(delta,   rank,0,m-rank,n);
                                        View sigma1(sigma[0],rank,0,m-rank,m);
                                        _BMD.mul(delta1,sigma1,serie[k]);                                        
                                        _BMD.mulin_right(Bperm, sigma[0]);
                                        for(size_t i=1;i<=std::min(k,max_degree);i++){
                                                View sigmak(sigma[i],rank,0,m-rank,m);
                                                _BMD.axpyin(delta1,sigmak,serie[k-i]);
                                                _BMD.mulin_right(Bperm, sigma[i]);
                                        }
                                        _BMD.mulin_right(Bperm, delta);
                                }
                                //std::cout<<"******** k="<<k<<std::endl;
#ifdef __DEBUG_MBASIS
                                std::cout<<sigma<<std::endl;
                                delta.write(std::cout,Tag::FileFormat::Maple);
                                std::cout<<std::endl;                                
#endif
                                BlasMatrix<Field> delta_copy(delta);
                                //delta_copy.write(std::cout,Tag::FileFormat::Plain);
                                
#define NEWELIM
#ifdef NEWELIM
                                // Compute a column reduced basis of the nullspace of delta
                                // -> [ L2 I ]                                 
                                rank= FFPACK::PLUQ (field(), FFLAS::FflasNonUnit, m, n,
                                                    delta_copy.getWritePointer(),delta_copy.getStride(),
                                                    Qt.getWritePointer(), P.getWritePointer());
#ifdef __DEBUG_MBASIS
                                delta_copy.write(std::cout,Tag::FileFormat::Maple);
                                std::cout<<std::endl;                                
#endif
                                View L1(delta_copy,0   ,0,  rank,rank);
                                View L2(delta_copy,rank,0,m-rank,rank);
                                FFLAS::ftrsm(field(),FFLAS::FflasRight,FFLAS::FflasLower,
                                             FFLAS::FflasNoTrans,FFLAS::FflasUnit,
                                             m-rank,rank, field().mOne, L1.getPointer(),L1.getStride(), L2.getWritePointer(),L2.getStride());
#ifdef __DEBUG_MBASIS
                                delta_copy.write(std::cout,Tag::FileFormat::Maple);
                                std::cout<<std::endl;                                
#endif

                                
#else
                                LQUPMatrix<Field> LQUP(delta_copy,P,Qt);
                                // Get L from LQUP
                                TriangularBlasMatrix<Field> L(field(), m, m, Tag::Shape::Lower, Tag::Diag::Unit);
                                LQUP.getL(L);
                                rank=LQUP.getRank();  // the first rank entries of Qt give the pivot row
                                // inverse L in-place (COULD BE IMPROVED -> only compute the left kernel by trsm)
                                //FFPACK::ftrtri(field(),FFLAS::FflasLower,FFLAS::FflasUnit,m,L.getPointer(),L.getStride());
                                View L1(L,0   ,0,  rank,rank);
                                View L2(L,rank,0,m-rank,rank);
                                FFLAS::ftrsm(field(),FFLAS::FflasRight,FFLAS::FflasLower,
                                             FFLAS::FflasNoTrans,FFLAS::FflasUnit,
                                             m-rank,rank, field().mOne, L1.getPointer(),m, L2.getWritePointer(),m);
#endif
                                
                                // update sigma by L^(-1) (rank sensitive -> use only the left kernel basis)
                                for(size_t i=0;i<=std::min(k,max_degree);i++){
                                        // NEED TO APPLY Qt to sigma[i]
                                        _BMD.mulin_right(Qt, sigma[i]);                                        
                                        View S1(sigma[i],0,0,rank,m);
                                        View S2(sigma[i],rank,0,m-rank,m);
                                        _BMD.axpyin(S2,L2,S1);
                                        //_BMD.mulin_right(L,sigma[i]);
                                }
#ifdef __DEBUG_MBASIS
                                std::cout<<"Qt=";
                                Qt.write(std::cout,false);
                                std::cout<<"\nP=";
                                P.write(std::cout,false);

                                std::cout<<std::endl<<"rank="<<rank<<" Qt size: "<<Qt.getSize()<<std::endl;
#endif
#if 0
                                size_t dmax=0, smax=0;
                                // update: the row-degree, the shifted row-degree,
                                //         the max pivot degree and the maximum row degree
                                for (size_t i=0;i<rank;++i) {
                                        //std::cout<<"Qt["<<i<<"]="<<Qt.getPointer()[i]<<std::endl;
                                        dmax=std::max(dmax, degree[Qt[i]]);
                                        smax=std::max(smax, shift [Qt[i]]);                                        
                                        degree[Qt[i]]++;
                                        shift [Qt[i]]++;
                                }
                                //std::cout<<"dmax:"<<dmax<<std::endl;
                                max_degree=std::max(max_degree,dmax+1);
                                //std::cout<<"max degree:"<<max_degree<<std::endl;
                                for (size_t i=rank;i<m;i++){
                                        degree[Qt[i]]=std::max(dmax, degree[Qt[i]]);
                                        //THIS LINE IS NOT NEEDED -> shift [Qt[i]]=max(smax, shift [Qt[i]]);
                                }
#else

                                size_t dmax=0;
                                Givaro::ZRing<size_t> Zint;
                                BlasMatrixDomain<Givaro::ZRing<size_t> > TTT(Zint);

                                TTT.mulin_right(Qt, degree);
                                TTT.mulin_right(Qt, shift);
#ifdef __DEBUG_MBASIS
                                std::cout<<"Qt.Bp.shift=";std::copy(shift.begin(),shift.end(),std::ostream_iterator<size_t>(std::cout,","));std::cout<<std::endl;
                                std::cout<<"Qt.Bp.degree=";std::copy(degree.begin(),degree.end(),std::ostream_iterator<size_t>(std::cout,","));std::cout<<std::endl;                                

                                std::cout<<std::endl;
#endif
#endif
                                
                                for (size_t i=0;i<rank;++i) {
                                        dmax=std::max(dmax, degree[i]);
                                        degree[i]++;
                                        shift [i]++;
                                }
                                for (size_t i=rank;i<m;i++){
                                        degree[i]=std::max(degree[i],dmax);
                                }
                                max_degree=std::min(order,std::max(max_degree,dmax+1));                     
#ifdef __DEBUG_MBASIS
                                std::cout<<"max degree:"<<max_degree<<std::endl;
                                std::cout<<"Qt.Bp.shift=";std::copy(shift.begin(),shift.end(),std::ostream_iterator<size_t>(std::cout,","));std::cout<<std::endl;
                                std::cout<<"Qt.Bp.degree=";std::copy(degree.begin(),degree.end(),std::ostream_iterator<size_t>(std::cout,","));std::cout<<std::endl;                                
                                
#endif
                                // shift the pivot row of sigma by x
                                for (int l=max_degree-1;l>=0;l--)
                                        for (size_t i=0;i<rank;i++)
                                                for (size_t j=0;j<m;j++)
                                                        sigma.ref(i,j,l+1)=sigma.ref(i,j,l);
                                                        //sigma.ref(Qt[i],j,l+1)=sigma.ref(Qt[i],j,l);
                                for (size_t i=0;i<rank;i++)
                                        for (size_t j=0;j<m;j++)
                                                sigma.ref(i,j,0)=field().zero;
#ifdef __DEBUG_MBASIS
                                std::cout<<"max degree="<<max_degree<<std::endl<<std::endl;
                                std::cout<<"F"<<k<<":="<<sigma<<std::endl<<"******************"<<std::endl;
                                std::cout<<std::endl;
#endif
                                // update Early Termination
                                //_EarlyStop.update(rank,shift); 
#ifdef __CHECK_MBASIS
                                std::cout<<"MBASIS: order "<<k<<check_orderbasis(field(),sigma,serie,k+1,0,shift)<<std::endl;
#endif
                                
                                _EarlyStop.update(m-rank,shift); // codimension (m-rank) seems better
                        }
                        
                        if (_EarlyStop.terminated()) { 
                                std::cout<<"OrderBasis: Early Termination at :"<<k<<"/"<<order<<std::endl;
                        }
                         
                        sigma.resize(max_degree+1);
                        return max_degree;
                }



                inline size_t twoValuation(size_t x){
                        size_t i=0;
                        while (x!=0 && !(x&0x1)){
                                i++;
                                x>>=1;
                        }
                        return i;
                }

                template<class Polynomial1>
                void update_sigma(size_t m, size_t n, std::list<Polynomial1*>& L, size_t k){
                        Polynomial1 *P1,*P2,*P3;
                        for(size_t i=0;i<k;i++){
                                P2=L.back();L.pop_back();
                                P3=L.back();L.pop_back();
                                P1 = new Polynomial1(field(),m,n,P2->size()+P3->size()-1);
                                _PMD.mul(*P1,*P2,*P3);
                                L.push_back(P1);
                                delete P2; delete P3;
                        }
                }

                // serie must have exactly order elements (i.e. its degree = order-1)
                // sigma can have at most order+1 elements (i.e. its degree = order)
                // Algorithm from [Giorgi, Lebreton ISSAC'2014]
                template<typename PMatrix1, typename PMatrix2>
                void oPM_Basis(PMatrix1             &sigma,
                               const PMatrix2       &serie,
                               size_t                order,
                               std::vector<size_t>       &shift)
                {
                        size_t m,n,k,l,lp;
                        m=sigma.rowdim();
                        n=sigma.coldim();
                        k=serie.coldim();

                        // log of the order
                        size_t log_order=integer(order).bitsize();

                        //  leaf size of the recursive PM_Basis algorithm (must be a power of 2)
                        size_t log_ord = MBASIS_THRESHOLD_LOG;
                        size_t ord     = std::min(size_t(1)<<log_ord ,order);

                        // prepare the storage for each serie update
                        std::vector<PMatrix2*> L_serie(log_order+1);
                        for (size_t i=log_ord;i<log_order;i++)
                                L_serie[i]= new PMatrix2(field(),n,k,(1<<i));
                        L_serie[log_order]=const_cast<PMatrix2*>(&serie);

                        typedef typename PMatrix2::const_view cview;
                        std::list<PMatrix1*>  L_sigma;
                        PMatrix1*         sigmak;
                        cview             seriek;
                        typedef HalflineMPDomain<Field,PMatrix2,PMatrix1,PMatrix2> HFMPD;
                        std::list<HFMPD*> L_mp;
                        typename std::list<HFMPD*>::iterator iter, t_iter;

                        // // Reset Early Termination
                        _EarlyStop.reset();

                        sigmak = new PMatrix1(field(),m,n,ord+1);
                        seriek = serie.at(0,ord-1);
                        M_Basis(*sigmak, seriek, ord, shift);
                        L_sigma.push_back(sigmak);
                        size_t sss=0;
                        for(size_t k=ord;k<order &&  !_EarlyStop.terminated();k+=ord,sss++){
                                //std::cout<<"------------ order="<<k<<std::endl;
                                l  = twoValuation(k);
                                lp = twoValuation(k-(1<<l)); lp=(lp==0?log_order:lp);
                                // compute next element in the original serie and
                                // update all subsequent computed series
                                for(iter=L_mp.begin(); iter!=L_mp.end(); ){
                                        (*iter)->update(ord);
                                        t_iter=iter;
                                        ++iter;

                                        if ((*t_iter)->terminated()){
                                                delete *t_iter;
                                                L_mp.erase(t_iter);
                                        }
                                }

                                // compute the serie update
                                //seriek = const_cast<const PMatrix2*>(L_serie[lp])->at(0,min(order,(1ULL<<(l+1)))-1);

                                //size_t update_max=min(order,1ULL<<(l+1));
                                //_PMD.midproductgen(*L_serie[l],*L_sigma.back(), seriek, true, (1ULL<<l)+1,update_max);

 
                                /*
                                  cout<<"---------------"<<endl;
                                  cout<<"MP "<<(1ULL<<l)<<"x"<<(1ULL<<(l+1))<<endl;
                                  cout<<*(L_sigma.back())<<endl;
                                  cout<<lp<<" ----- "<<endl<<*(L_serie[lp])<<endl;
                                  cout<<"---------------"<<endl;
                                */

                                //L_mp.push_back(new HFMPD(field(),*(L_serie[l]),*(L_sigma.back()), seriek, 1ULL<<l));
                                L_mp.push_back(new HFMPD(field(),*(L_serie[l]),*(L_sigma.back()),*(L_serie[lp]), 1ULL<<l));

                                // compute the new sigma
                                sigmak = new PMatrix(field(),m,n,ord+1);
                                size_t step= std::min(ord,order-k); // needed if order%ord <> 0
                                seriek = const_cast<const PMatrix2*>(L_serie[l])->at(0,step-1);

                                // compute the next "step" elements of the serie L_serie[l]
                                L_mp.back()->update(ord);

                                M_Basis(*sigmak, seriek, step, shift);
                                L_sigma.push_back(sigmak);

                                // update the sigma list
                                update_sigma(m, n, L_sigma, twoValuation(k/ord+1));
                        }
                        // get the product of all sigma
                        if (L_sigma.size()>1){
                                update_sigma(m, n, L_sigma, L_sigma.size()-2);
                                PMatrix1 *s1,*s2;
                                s1= L_sigma.back();L_sigma.pop_back();
                                s2= L_sigma.back();L_sigma.pop_back();
                                _PMD.mul(sigma,*s1,*s2);
                                delete s1; delete s2;
                        }
                        else {
                                sigma.resize(L_sigma.back()->size());
                                sigma.copy(*L_sigma.back(),0,L_sigma.back()->size()-1);
                                delete L_sigma.back();
                        }
                        for (size_t i=0;i<log_order;i++)
                                delete L_serie[i];

                        // Info about early termination
                        //if (_EarlyStop.terminated())
                        //        cout<<"Early termination at order "<<sss<<" ("<<order<<")"<<endl;
                }

#ifdef LOW_MEMORY_PMBASIS
                // serie must have exactly order elements (i.e. its degree = order-1)
                // sigma can have at most order+1 elements (i.e. its degree = order)
                // !!! sigma is not allocated apriori !!!
                template<typename PMatrix1, typename PMatrix2>
                size_t PM_Basis_low(PMatrix1*                &sigma_ptr,
                                    const PMatrix2           *serie_ptr,
                                    size_t                    order,
                                    std::vector<size_t>       &shift)
                {

#ifdef PROFILE_PMBASIS
                        //std::cout<<"Start PM-Basis : "<<order<<" ("<<_idx<<"/"<<_target<<")] : "<<std::endl;//MEMINFO2<<std::endl;
                        if (_target==0) _target=order;
                        if (!_started) {_started=true; _start = std::chrono::system_clock::now();}
                        std::chrono::time_point<std::chrono::system_clock> _chrono_start=std::chrono::system_clock::now();
#endif
                        
                        if (order <= MBASIS_THRESHOLD) {
#if defined (PROFILE_PMBASIS) or defined(__CHECK_PMBASIS)
                                _idx+=order;
#endif
                                sigma_ptr = new PMatrix1(field(),serie_ptr->rowdim(),serie_ptr->rowdim(),order+1);
                                size_t res= M_Basis(*sigma_ptr, *serie_ptr, order, shift);
                                delete serie_ptr;
                                return res;
                        }
                        else {
#ifdef PROFILE_PMBASIS
                                Timer chrono;
                                chrono.start();
#endif
                                size_t ord1,ord2,d1,d2;
                                //ord1 = order>>1;
                                //ord2 = order-ord1; // ord1+ord2=order
                                ord2 = order>>1;
                                ord1 = order-ord2; // ord1+ord2=order
                                size_t m,n,k;
                                m=serie_ptr->rowdim();
                                n=serie_ptr->rowdim();
                                k=serie_ptr->coldim();
                                integer p;

                                // first recursive call
                                PMatrix1 *sigma1_ptr, *sigma2_ptr;
                                PMatrix2 *serie1_ptr, *serie2_ptr;

                                // Allocate serie1
                                serie1_ptr= new PMatrix2(field(),n,k,ord1);                                
#ifdef MEM_PMBASIS
                                std::cerr<<"[PM-Basis ("<<order<<") "<<_idx<<"/"<<_target<<"] [ALLOC Serie1] -> "<<MB(serie1_ptr->realmeminfo())<<"Mo"<<MEMINFO2<<std::endl;
#endif
                                serie1_ptr->copy(*serie_ptr,0,ord1-1);
                                d1 = PM_Basis_low(sigma1_ptr, serie1_ptr, ord1, shift);
                                // no more needed
                                // delete serie1_ptr; 
                                
#ifdef MEM_PMBASIS
                                std::cerr<<"[PM-Basis ("<<order<<") "<<_idx<<"/"<<_target<<"] [DEL Serie1] -> "<<MEMINFO2<<std::endl;
#endif


                                if (_EarlyStop.terminated()){
                                        sigma_ptr=sigma1_ptr;
                                        delete serie_ptr;
                                        return d1;
                                }

                                // Allocate serie2
                                serie2_ptr=new PMatrix2(field(),n,k,ord2);//serie2 size=ord1+1 -> midproduct)
#ifdef MEM_PMBASIS
                                std::cerr<<"[PM-Basis ("<<order<<") "<<_idx<<"/"<<_target<<"] [ALLOC Serie2] -> "<<MB(serie2_ptr->realmeminfo())<<"Mo"<<MEMINFO2<<std::endl;
#endif
                                
                                _PMD.midproductgen(*serie2_ptr, *sigma1_ptr, *serie_ptr, true, ord1+1,ord1+ord2);
#ifndef __CHECK_PMBASIS
                                delete serie_ptr; // the initial serie is no more needed (except with checking pmbasis)
#endif         
                                // second recursive call                                                                
                                d2 = PM_Basis_low(sigma2_ptr, serie2_ptr, ord2, shift);
                                // no more needed
                                // delete serie2_ptr;
#ifdef MEM_PMBASIS
                                std::cerr<<"[PM-Basis ("<<order<<") "<<_idx<<"/"<<_target<<"] [DEL Serie2] -> "<<MEMINFO2<<std::endl;
#endif                                
                                // compute the result
                                sigma_ptr = new PMatrix1(field(),m,n,d1+d2+1);
                                //sigma_ptr = new PMatrix1(field(),m,n,order+1);                                
#ifdef MEM_PMBASIS
                                std::cerr<<"[PM-Basis ("<<order<<") "<<_idx<<"/"<<_target<<"] [ALLOC Sigma] -> "<<MB(sigma_ptr->realmeminfo())<<"Mo"<<MEMINFO2<<std::endl;
#endif                                
                                _PMD.mul(*sigma_ptr, *sigma2_ptr, *sigma1_ptr, d1+d2);
                                //sigma_ptr->resize(d1+d2+1);                                
                                delete sigma1_ptr;
                                delete sigma2_ptr;
#ifdef MEM_PMBASIS
                                std::cerr<<"[PM-Basis ("<<order<<") "<<_idx<<"/"<<_target<<"] [DEL Sigma 1/2] -> "<<MEMINFO2<<std::endl;
#endif

                                
#ifdef PROFILE_PMBASIS
                                //chrono.stop();
                                //std::cout<<"      -> basis product "<<sigma1.size()<<"x"<<sigma2.size()<<" = "<<d1+d2+1<<" -->"<<chrono.usertime()<<MEMINFO2<<std::endl;
#endif

#ifdef __CHECK_PMBASIS
                                if(order >= __CHECK_PMBASIS_THRESHOLD){
                                        std::cerr<<"PMBASIS: order "<<order<<check_orderbasis(field(),*sigma_ptr,*serie_ptr,order,((_idx/order)+1)&1,shift)<<std::endl;                                        
                                }
                                delete serie_ptr;
#endif
#ifdef PROFILE_PMBASIS
                                chrono.stop();
                                _end = std::chrono::system_clock::now();                                
                                std::chrono::duration<double> elapsed_beginning = _end-_start;
                                std::chrono::duration<double> elapsed_comp      = _end-_chrono_start;

                                double magicnumber=double(_target)/double(order)*log(double(_target)/double(order))/log(2.);
                                double tcomp = elapsed_comp.count();
                                double telap = elapsed_beginning.count();
                                
                                _eta=(_eta!=0.0?std::min(_eta,tcomp*magicnumber):tcomp*magicnumber);
                                std::cerr<<"[PM-Basis : "<<order<<" ("<<_idx<<"/"<<_target<<")] : "<<chrono.usertime()
                                         << " (ETA: "<< telap<<"s / "<<_eta<<"s)"<<MEMINFO2<<std::endl;
                                chrono.clear();chrono.start();
#endif


                                return d1+d2;
                        }
                }
#endif // LOW_MEMORY_PMBASIS

    // new mbasis: output is in s-ordered weak Popov form
    std::vector<size_t> mbasis(
        PMatrix &approx,
        const PMatrix &series,
        const size_t order,
        const std::vector<int> &shift=std::vector<int>(),
        bool resUpdate=false );
    
    // new pmbasis: output is in s-ordered weak Popov form
    std::vector<size_t> pmbasis(
        PMatrix &approx,
        const PMatrix &series,
        const size_t order,
        const std::vector<int> &shift=std::vector<int>(),
	const size_t threshold=32 );
        
    // pmbasis with output in s-Popov form
    std::vector<size_t> popov_pmbasis(
        PMatrix &approx,
        const PMatrix &series,
        const size_t order,
        const std::vector<int> &shift=std::vector<int>(),
	const size_t threshold=32 );

        }; // end of class OrderBasis

        
        typedef Givaro::Modular<RecInt::ruint128,RecInt::ruint256>   MYRECINT;
        template<>
		size_t OrderBasis<MYRECINT,EarlyTerm<(size_t) -1> >::M_Basis(PolynomialMatrix<PMType::polfirst,PMStorage::plain, MYRECINT>            &sigma,
                                                            const PolynomialMatrix<PMType::polfirst,PMStorage::plain, MYRECINT>      &serie,
                                                            size_t                 order,
                                                            std::vector<size_t>   &shift)
        {
                Givaro::Integer p; field().cardinality(p);
                typedef Givaro::Modular<Givaro::Integer> NewField;
                NewField F(p);
                OrderBasis<NewField > SB(F);
                typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain, NewField> NewMatrix;
                
                NewMatrix sigma1(F,sigma.rowdim(),sigma.coldim(),order+1);
                NewMatrix serie1(F,serie.rowdim(),serie.coldim(),order);
                serie1.copy(serie,0,order-1);

                //std::cout<<"Serie: "<<serie<<std::endl;
                //std::cout<<"Serie1: "<<serie1<<std::endl;

                size_t d= SB.M_Basis(sigma1,serie1,order,shift);
                sigma.copy(sigma1,0,d);
                
                //std::cout<<"Sigma1: "<<sigma1<<std::endl;
                //std::cout<<"Sigma: "<<sigma<<std::endl;


                return d;
        }
        
/* TODO: vneiger: check if correct for both resUpdate==false and
 * resUpdate==true; a threshold should be defined here to know where to switch
 * between continuous update of the residual and full recomputation... should
 * maybe even switch within the same calls, depending on the current reached
 * order?  */
/** Algorithm M-Basis-One as detailed in Section 3 of
 *  [Jeannerod, Neiger, Villard. Fast Computation of approximant bases in
 *  canonical form. Preprint, 2017]
 **/
/** Input:
 *   - approx: m x m square polynomial matrix, approximation basis
 *   - series: m x n polynomial matrix of size = order, series to approximate
 *   - order: positive integer, order of approximation
 *   - shift: m-tuple of degree shifts (acting on columns of approx)
 *   - resUpdate: FIXME
 **/
/** Action:
 *   - Compute and store in 'approx' an shift-ordered weak Popov
 *   approximation basis for (series,order)
 *  Note: if order=1, then approx is in shift-Popov form
 **/
/** Output: shifted minimal degree of (series,order),
 * which is equal to the diagonal degrees of approx **/
/** Complexity: O(m^w order^2) **/
template<class Field, class ET>
std::vector<size_t> OrderBasis<Field,ET>::mbasis( OrderBasis<Field,ET>::PMatrix &approx, const OrderBasis<Field,ET>::PMatrix &series, const size_t order, const std::vector<int> &shift, bool resUpdate )
{
    const size_t m = series.rowdim();
    const size_t n = series.coldim();
    typedef BlasSubmatrix<typename OrderBasis<Field,ET>::MatrixP::Matrix> View;

    // initialize approx to the identity matrix
    approx.resize(0); // to put zeroes everywhere.. FIXME may be a better way to do it but it seems approx.clear() fails
    size_t appsz = 1;
    approx.resize(appsz);
    for ( size_t i=0; i<m; ++i )
        approx.ref(i,i,0) = 1;

    // initial shifted row degrees = shift
    std::vector<int> rdeg( shift );
    // initial shifted minimal degree = (0,...,0)
    std::vector<size_t> mindeg( m, 0 );

    // set residual to input series
    OrderBasis<Field,ET>::PMatrix res( this->field(), m, n, series.size() );
    res.copy( series );

    for ( size_t ord=0; ord<order; ++ord )
    {
        //At the beginning of iteration 'ord',
        //   - approx is an order basis, shift-ordered weak Popov,
        //   for series at order 'ord'
        //   - the shift-min(shift) row degrees of approx are rdeg.
        //   - the max degree in approx is <= appsz

        // Here we follow [Algorithm M-Basis-1] in the above reference

        // coefficient of degree 'ord' of residual, which we aim at cancelling
        typename OrderBasis<Field,ET>::MatrixP::Matrix res_const( approx.field(), m, n );
        if ( resUpdate ) // res_const is coeff of res of degree ord
            res_const = res[ord];
        else // res_const is coeff of approx*res of degree ord
        {
            for ( size_t d=0; d<appsz; ++d )
                this->_BMD.axpyin( res_const, approx[d], res[ord-d] ); // note that d <= appsz-1 <= ord
        }

        // permutation for the stable sort of the shifted row degrees
        std::vector<size_t> perm_rdeg( m );
        iota( perm_rdeg.begin(), perm_rdeg.end(), 0 );
        stable_sort(perm_rdeg.begin(), perm_rdeg.end(),
                [&](const size_t& a, const size_t& b)->bool
                {
                    return (rdeg[a] < rdeg[b]);
                } );

        // permute rows of res_const accordingly
        std::vector<size_t> lperm_rdeg( m ); // LAPACK-style permutation
        FFPACK::MathPerm2LAPACKPerm( lperm_rdeg.data(), perm_rdeg.data(), m );
        BlasPermutation<size_t> pmat_rdeg( lperm_rdeg );
        this->_BMD.mulin_right( pmat_rdeg, res_const );

        // compute PLUQ decomposition of res_const
        BlasPermutation<size_t> P(m), Q(n);
        size_t rank = FFPACK::PLUQ( res_const.field(), FFLAS::FflasNonUnit, //FIXME TODO investigate see below ftrsm
                        m, n, res_const.getWritePointer(), res_const.getStride(),
                        P.getWritePointer(), Q.getWritePointer() );

        // compute a part of the left kernel basis of res_const:
        // -Lbot Ltop^-1 , stored in Lbot
        // Note: the full kernel basis is [ -Lbot Ltop^-1 | I ] P
        View Ltop( res_const, 0, 0, rank, rank ); // top part of lower triangular matrix in PLUQ
        View Lbot( res_const, rank, 0, m-rank, rank ); // bottom part of lower triangular matrix in PLUQ
        FFLAS::ftrsm( approx.field(), FFLAS::FflasRight, FFLAS::FflasLower,
                  FFLAS::FflasNoTrans, FFLAS::FflasUnit, // FIXME TODO works only if nonunit in PLUQ and unit here; or converse. But not if consistent...?????? investigate
                  m-rank, rank, approx.field().mOne,
                  Ltop.getPointer(), Ltop.getStride(),
                  Lbot.getWritePointer(), Lbot.getStride() );

        // Prop: this "kernel portion" is now stored in Lbot.
        //Then const_app = perm^{-1} P^{-1} [ [ X Id | 0 ] , [ Lbot | Id ] ] P perm
        //is an order basis in rdeg-Popov form for const_res at order 1
        // --> by transitivity,  const_app*approx will be a shift-ordered
        // weak Popov approximant basis for (series,ord+1)

        // A. update approx basis, first steps:
        for ( size_t d=0; d<appsz; ++d )
        {
            // 1. permute rows: approx = P * perm * approx
            this->_BMD.mulin_right( pmat_rdeg, approx[d] );
            this->_BMD.mulin_right( P, approx[d] );
            // 2. multiply by constant: appbot += Lbot apptop
            View apptop( approx[d], 0, 0, rank, m );
            View appbot( approx[d], rank, 0, m-rank, m );
            this->_BMD.axpyin( appbot, Lbot, apptop );
        }

        // permute row degrees accordingly
        std::vector<size_t> lperm_p( P.getStorage() ); // Lapack-style permutation P
        std::vector<size_t> perm_p( m ); // math-style permutation P
        FFPACK::LAPACKPerm2MathPerm( perm_p.data(), lperm_p.data(), m ); // convert

        // B. update shifted row degree, shifted minimal degree,
        // and new approximant basis size using property: deg(approx) = max(mindeg)
        for ( size_t i=0; i<rank; ++i )
        {
            ++rdeg[perm_rdeg[perm_p[i]]];
            ++mindeg[perm_rdeg[perm_p[i]]];
        }
        appsz = 1 + *max_element(mindeg.begin(),mindeg.end());
        approx.resize( appsz );

        // A. update approx basis:
        // 3. multiply first rank rows by X...
        for ( size_t d=appsz-1; d>0; --d )
            for ( size_t i=0; i<rank; ++i )
                for ( size_t j=0; j<m; ++j )
                    approx.ref(i,j,d) = approx.ref(i,j,d-1);
        // 4. ... and approx[0]: first rank rows are zero
        for ( size_t i=0; i<rank; ++i )
            for ( size_t j=0; j<m; ++j )
                approx.ref(i,j,0) = 0;
        // 5. permute the rows again: approx = perm^{-1} * P^{-1} * approx
        P.Invert();
        pmat_rdeg.Invert();
        for ( size_t d=0; d<appsz; ++d )
        {
            this->_BMD.mulin_right( P, approx[d] );
            this->_BMD.mulin_right( pmat_rdeg, approx[d] );
        }

        // TODO work on resUpdate = True; now assuming False
        if ( resUpdate )
        {
            //if ( resUpdate )
            //{
            //    for ( size_t d=ord; d<res.size(); ++d )
            //        this->_BMD.mulin_right( pmat_rdeg, res[d] );
            //}
            // update residual: do same operations as on approx
            // update residual: 1/ permute all the rows; multiply by constant
            // note: to simplify later multiplication by X, we permute rows of res[ord]
            //(but we don't compute the zeroes in the other rows)
            this->_BMD.mulin_right( P, res[ord] ); // permute rows by P
            for ( size_t d=ord+1; d<res.size(); ++d )
                this->_BMD.mulin_right( P, res[d] ); // permute rows by P

            for ( size_t d=ord+1; d<res.size(); ++d )
            {
                // multiply by constant: resbot += Lbot restop
                View restop( res[d], 0, 0, rank, n );
                View resbot( res[d], rank, 0, m-rank, n );
                this->_BMD.axpyin( resbot, Lbot, restop );
            }

            // update residual: 2/ multiply first rank rows by X...
            for ( size_t d=res.size()-1; d>ord; --d )
                for ( size_t i=0; i<rank; ++i )
                    for ( size_t j=0; j<n; ++j )
                        res.ref(i,j,d) = res.ref(i,j,d-1);
        }
    }
    return mindeg;
}

/* TODO: vneiger: needs better handling of the threshold; note that it might
 * differ for different fields (how to do this properly? threshold as attribute
 * of the class?) */
/** Algorithm PM-Basis as detailed in Section 2.2 of
 *  [Giorgi, Jeannerod, Villard. On the Complexity 
 *  of Polynomial Matrix Computations. ISSAC 2003]
 **/
/** Input:
 *   - approx: m x m square polynomial matrix, approximation basis
 *   - series: m x n polynomial matrix of degree < order, series to approximate
 *   - order: positive integer, order of approximation
 *   - shift: degree shift on the cols of approx
 *   - threshold: depth for leaves of recursion (when the current order reaches threshold, apply mbasis)
 **/
/** Action:
 *   - Compute and store in 'approx' a shifted ordered weak Popov
 *   approximation basis for (series,order,shift)
 **/
/* Return:
 * - the shifted minimal degree for (series,order)
 */
/** Output: shifted row degrees of the computed approx **/
/** Complexity: O(m^w M(order) log(order) ) **/
template<class Field, class ET>
std::vector<size_t> OrderBasis<Field,ET>::pmbasis(
    OrderBasis<Field,ET>::PMatrix &approx,
    const OrderBasis<Field,ET>::PMatrix &series,
    const size_t order,
    const std::vector<int> &shift,
    const size_t threshold )
{
    if ( order <= threshold )
    {
        std::vector<size_t> mindeg = mbasis( approx, series, order, shift );
        return mindeg;
    }
    else
    {
        size_t m = series.rowdim();
        size_t n = series.coldim();
        size_t order1,order2;
        order1 = order>>1; // order1 ~ order/2
        order2 = order - order1; // order2 ~ order/2, order1 + order2 = order
        std::vector<size_t> mindeg( m );

        OrderBasis<Field,ET>::PMatrix approx1( this->field(), m, m, 0 );
        OrderBasis<Field,ET>::PMatrix approx2( this->field(), m, m, 0 );

        {
            OrderBasis<Field,ET>::PMatrix res1( this->field(), m, n, order1 ); // first residual: series truncated mod X^order1
            res1.copy( series, 0, order1-1 );
            mindeg = pmbasis( approx1, res1, order1, shift ); // first recursive call
        } // end of scope: res1 is deallocated here
        {
            std::vector<int> rdeg( shift ); // shifted row degrees = mindeg + shift
            for ( size_t i=0; i<m; ++i )
                rdeg[i] += mindeg[i];
            OrderBasis<Field,ET>::PMatrix res2( series.field(), m, n, order2 ); // second residual: midproduct 
            this->_PMD.midproductgen( res2, approx1, series, true, order1+1, order1+order2 ); // res2 = (approx1*series / X^order1) mod X^order2
            std::vector<size_t> mindeg2( m );
            mindeg2 = pmbasis( approx2, res2, order2, rdeg ); // second recursive call
            for ( size_t i=0; i<m; ++i )
                mindeg[i] += mindeg2[i];
        } // end of scope: res2 is deallocated here
        
        // for PMD.mul we need the size to be the sum (even though we have a better bound on the output degree)
        //approx.resize( approx1.size()+approx2.size()-1 );
        // in fact, deg(approx) = max(mindeg)  (which is indeed the sum of the mindegs for approx1 and approx2)
        approx.resize( 1 + *max_element( mindeg.begin(), mindeg.end() ) );
        this->_PMD.mul( approx, approx2, approx1 );
        return mindeg;
    }

}

/* TODO: vneiger: introduce correct improvement in "uniform" case (what case
 * exactly? should be "column degree = mindeg"? or is the one mentioned below
 * correct?) */
/** Algorithm Popov-PM-Basis as detailed in
 *  [Jeannerod, Neiger, Villard. Fast Computation of approximant bases in
 *  canonical form. Preprint, 2017]
 **/
/** Input:
 *   - approx: m x m square polynomial matrix, approximation basis
 *   - series: m x n polynomial matrix of degree < order, series to approximate
 *   - order: positive integer, order of approximation
 *   - shift: degree shift on the cols of approx
 *   - threshold: depth for leaves of recursion (when the current order reaches threshold, apply mbasis)
 **/
/** Action:
 *   - Compute and store in 'approx' the shift-Popov approximation basis for (series,order)
 **/
/** Output: shifted minimal degree for (series,order) **/
/** Complexity: roughly O(m^w M(order) log(order) ) **/
template<class Field, class ET>
std::vector<size_t> OrderBasis<Field,ET>::popov_pmbasis(
    OrderBasis<Field,ET>::PMatrix &approx,
    const OrderBasis<Field,ET>::PMatrix &series,
    const size_t order,
    const std::vector<int> &shift,
    const size_t threshold )
{
    // 1. compute shift-ordered weak Popov approximant basis
    size_t m = series.rowdim();
    std::vector<size_t> mindeg = pmbasis( approx, series, order, shift, threshold );
    
    // 2. compute -mindeg-ordered weak Popov approximant basis

    //// Note: if mindeg+shift is uniform, (????)
    //// then approx is almost in -mindeg-??? form
    //bool uniform=true;
    //for ( size_t i=0; i<m-1; ++i )
    //{
    //    if ( (int)mindeg[i] + shift[i] != (int)mindeg[i+1] + shift[i+1]) {
    //        uniform = false;
    //    }
    //}

    //if ( !uniform )
    //{
        std::vector<int> mindegshift( m );
        for ( size_t i=0; i<m; ++i )
            mindegshift[i] = - (int)mindeg[i];
        pmbasis( approx, series, order, mindegshift, threshold );
    //}

    // 3. left-multiply by inverse of -mindeg-row leading matrix
    // Note: cdeg(approx) = mindeg
    typename OrderBasis<Field,ET>::MatrixP::Matrix lmat( this->field(), m, m );
    for ( size_t i=0; i<m; ++i )
    for ( size_t j=0; j<m; ++j )
        lmat.setEntry( i, j, approx.get( i, j, mindeg[j] ) );
    this->_BMD.invin( lmat );
    for ( size_t k=0; k<approx.size(); ++k ) {
        this->_BMD.mulin_right( lmat, approx[k] );
    }
    return mindeg;
}
        
} // end of namespace LinBox

// Local Variables:
// mode: C++ 
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
