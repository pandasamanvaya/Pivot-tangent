#include<seal/seal.h>
#include<iostream>
#include<vector>
#include<string>
#include<random>
#include<chrono>
#include<math.h>
#include "matplotlibcpp.h"

using namespace std;
using namespace seal;
namespace plt = matplotlibcpp;

class InvSqrt{
private:
	SEALContext* context = nullptr;

	/*Different Parameters for pivot-tangent method */

	// a = 1e-4; b = 1e3; d = 8; d_g = 7; err = 7e-3
	// double x_l = 1e-4, x_r = 1e3;
	// double k2 = 1.6876766119199442, k1 = 0.0855726219177246;
	// double x2 = 340.0359057573916, x1 = 0.13221763833040837;
	// double pivot = 0.3888343439376094;

	//a = 1e-4; b = 1e3; d = 9; d_g = 7; err = 1e-3
	// double x_l = 1e-4, x_r = 1e3;
	// double k2 = 1.6958489466846607, k1 = 0.070270068359375;
	// double x2 = 338.78048579006565, x1 = 0.0876277694155693;
	// double pivot = 0.2586564843344327;
	
	//a = 1e-4; b = 1e3; d = 9; d_g = 7; err = 1e-4
	double x_l = 1e-4, x_r = 1e3;
	double k2 = 1.6885913944245319, k1 = 0.08385114624023438;
	double x2 = 339.8933847174716, x1 = 0.12670372040311723;
	double pivot = 0.37277489383454;

	//a = 1e-4; b = 1e3; d = 7; d_g = 7; err = 7e-3
	// double x_l = 1e-4, x_r = 1e3;
	// double k2 = 1.6644849487774667, k1 = 0.1280319860839844;
	// double x2 = 343.66449980138054, x1 = 0.3111193884858657;
	// double pivot = 0.9052997579490315;

	//a = 1e-3; b = 750; d = 7; d_g = 5; err = 1e-3
	// double x_l = 1e-3, x_r = 750;
	// double k2 = 1.6480555806754493, k1 = 0.15725496246337894;
	// double x2 = 259.7224294229279, x1 = 0.3649351552907927;
	// double pivot = 1.053822602369088;

	//a = 1e-3; b = 750; d = 7; d_g = 5; err = 7e-3
	// double x_l = 1e-3, x_r = 750;
	// double k2 = 1.6286249562069153, k1 = 0.1909952459716797;
	// double x2 = 262.11481941471857, x1 = 0.5622392419520039;
	// double pivot = 1.6087584532823394;
	
	//a = 1e-3; b = 750; d = 8; d_g = 6; err = 1e-4
	// double x_l = 1e-3, x_r = 750;
	// double k2 = 1.6658932150882901, k1 = 0.12548308563232424;
	// double x2 = 257.58027649210305, x1 = 0.22345451540236524;
	// double pivot = 0.6506355565501228;

	double m1 = -0.5*k2*pow(x1, -1.5), c1 = 1.5*k2*pow(x1,-0.5);
	double m2 = -0.5*k2*pow(x2, -1.5), c2 = 1.5*k2*pow(x2,-0.5);

	Plaintext a, b, half, div_b, neg_p, err, M1, M2, C1, C2;
public:
	Encryptor* encryptor = nullptr;
	Decryptor* decryptor = nullptr;
	CKKSEncoder* encoder = nullptr;
	Evaluator* evaluator = nullptr;
	RelinKeys* relin_keys = nullptr;
	GaloisKeys* galois_keys = nullptr;
	double scale;
	size_t N; 
	InvSqrt(SEALContext& context, Encryptor& encryptor, Decryptor& decryptor, CKKSEncoder& encoder, Evaluator& evaluator,
				 double scale, RelinKeys& relin_keys)
	{
	    this->scale = scale;
	    this->encryptor = &encryptor;
	    this->decryptor = &decryptor;
	    this->encoder = &encoder;
	    this->evaluator = &evaluator;
	    this->relin_keys = &relin_keys;
	    // this->galois_keys = &galois_keys;
	    this->context = &context;

	    N = encoder.slot_count();
   	    this->encoder->encode(x_l, scale, a);
	    this->encoder->encode(x_r, scale, b);
	    this->encoder->encode(1/(x_r-x_l), scale, div_b);
	    this->encoder->encode(-pivot, scale, neg_p);
		this->encoder->encode(0.5, scale, half);
	    this->encoder->encode(m1, scale, M1);
		this->encoder->encode(m2, scale, M2);
		this->encoder->encode(c1, scale, C1);
		this->encoder->encode(c2, scale, C2);

	}

	double getBounds(int i=1)
	{
		if(i == 1)
			return k1;
		return k2;
	}

	void reEncrypt(Ciphertext& ct)
	{
		Plaintext temp;
		vector<double> v;
		decryptor->decrypt(ct, temp);
		encoder->decode(temp, v);
		encoder->encode(v, scale, temp);
		encryptor->encrypt(temp, ct);
		// cout << "depth = " << context->get_context_data(ct.parms_id())->chain_index() << "\n";
	}

	Ciphertext evalLine(Ciphertext x, Plaintext m, Plaintext c)
	{
		// cout << "line\n";
		evaluator->mod_switch_to_inplace(m, x.parms_id());
		evaluator->multiply_plain_inplace(x, m);
		evaluator->rescale_to_next_inplace(x);
		evaluator->mod_switch_to_inplace(c, x.parms_id());
		x.scale() = scale;
		evaluator->add_plain_inplace(x, c);
		return x;
	}

	void printVector(Ciphertext ct, int l=1, int index=-1)
	{
		Plaintext temp;
		vector<double> result;
		decryptor->decrypt(ct, temp);
		encoder->decode(temp, result);

		if(index != -1)
			// cout << result[index] << "\n";
			printf("%0.10f \n", result[index]);
		else{
			for(int i=0; i < l; i++)
				// cout << result[i] << " ";
				printf("%0.10f ", result[i]);
			cout << "\n";			
		}
	}

	int depth(Ciphertext c)
	{
		return context->get_context_data(c.parms_id())->chain_index();
	}

	Ciphertext signPolyEval(Ciphertext x, vector<Plaintext> coeff)
	{
		// cout << "Initial depth " << context->get_context_data(x.parms_id())->chain_index() << "\n";
	    // x^2
	    Ciphertext x_2;
	    evaluator->square(x, x_2);
	    evaluator->relinearize_inplace(x_2, *relin_keys);
	    evaluator->rescale_to_next_inplace(x_2);
	    // cout << "x: " << x.scale() << " x_2: " << x_2.scale() << "\n";
	    
	    // x^3
	    Ciphertext x_3;
	    evaluator->mod_switch_to_inplace(x, x_2.parms_id());
	    evaluator->multiply(x_2, x, x_3);
	    evaluator->relinearize_inplace(x_3, *relin_keys);
	    evaluator->rescale_to_next_inplace(x_3);
	    // cout << "x_2: " << x_2.scale() << " x_3: " << x_3.scale() << "\n";

	    //x^4
	    Ciphertext x_4;
	    evaluator->square(x_2, x_4);
	    evaluator->relinearize_inplace(x_4, *relin_keys);
	    evaluator->rescale_to_next_inplace(x_4);
	    // cout << "x_4: " << x_4.scale() << " x_2: " << x_2.scale() << "\n";

	    //x^5
	    Ciphertext x_5;
	    evaluator->mod_switch_to_inplace(x_2, x_3.parms_id());
	    evaluator->multiply(x_2, x_3, x_5);
	    evaluator->relinearize_inplace(x_5, *relin_keys);
	    evaluator->rescale_to_next_inplace(x_5);
	    // cout << "x_5: " << x_5.scale() << " x_3: " << x_3.scale() << "\n";

	    //x^7
	    Ciphertext x_7;
	    evaluator->multiply(x_3, x_4, x_7);
	    evaluator->relinearize_inplace(x_7, *relin_keys);
	    evaluator->rescale_to_next_inplace(x_7);
	    // cout << "x_7: " << x_7.scale() << "\n";

	    // cout << depth(x_7) <<"\n";
	    //Multiply constants
	    evaluator->mod_switch_to_inplace(coeff[1], x.parms_id());
	    evaluator->multiply_plain_inplace(x, coeff[1]);
	    evaluator->rescale_to_next_inplace(x);

	    evaluator->mod_switch_to_inplace(coeff[3], x_3.parms_id());
	    evaluator->multiply_plain_inplace(x_3, coeff[3]);
	    evaluator->rescale_to_next_inplace(x_3);

	    evaluator->mod_switch_to_inplace(coeff[5], x_5.parms_id());
	    evaluator->multiply_plain_inplace(x_5, coeff[5]);
	    evaluator->rescale_to_next_inplace(x_5);

	    evaluator->mod_switch_to_inplace(coeff[7], x_7.parms_id());
	    evaluator->multiply_plain_inplace(x_7, coeff[7]);
	    evaluator->rescale_to_next_inplace(x_7);
	    //c7*x^7 + c5*x^5 + c3*x^3 + c1*x
	    x.scale() = scale; x_3.scale() = scale; x_5.scale() = scale; x_7.scale() = scale;
	    evaluator->mod_switch_to_inplace(x, x_7.parms_id());
	    evaluator->mod_switch_to_inplace(x_3, x_7.parms_id());
	    evaluator->mod_switch_to_inplace(x_5, x_7.parms_id());
	    evaluator->add_inplace(x, x_3);
	    evaluator->add_inplace(x, x_5);
	    evaluator->add_inplace(x, x_7);

	    // cout << "Final depth " << context->get_context_data(x.parms_id())->chain_index() << "\n";

		return x;		
	}

	Ciphertext signApprox(Ciphertext x, int d_g, int d_f)
	{
		vector<double> df_coeff = {0, 35.0/16, 0, -35.0/16, 0, 21.0/16, 0, -5.0/16};
		vector<double> dg_coeff = {0, 4589.0/1024, 0, -16577.0/1024, 0, 25614.0/1024, 0, -12860.0/1024};

		vector<Plaintext> encode_df(8), encode_dg(8);
		for(int i=0; i < 8; i++){
			encoder->encode(df_coeff[i], scale, encode_df[i]);
			encoder->encode(dg_coeff[i], scale, encode_dg[i]);
		}
		for(int i=0; i < d_g; i++){
			if(depth(x) < 4)
				reEncrypt(x);
			x = signPolyEval(x, encode_dg);
			// x = reEncrypt(x);
			// printVector(x, 5);
		}
		for(int i=0; i < d_f; i++){
			if(depth(x) < 4)
				reEncrypt(x);
			x = signPolyEval(x, encode_df);
			// x = reEncrypt(x);
			// printVector(x, 5);
		}
		reEncrypt(x);
		return x;
	}

	Ciphertext comp(Ciphertext x, int d_g, int d_f)
	{
		Plaintext one;
		encoder->encode(1.0, scale, one);
		//X = (x-p)/b
		evaluator->mod_switch_to_inplace(neg_p, x.parms_id());
		evaluator->add_plain_inplace(x, neg_p);
		evaluator->mod_switch_to_inplace(div_b, x.parms_id());
		evaluator->multiply_plain_inplace(x, div_b);
		evaluator->rescale_to_next_inplace(x);
		// cout << "normalise\n";
		x = signApprox(x, d_g, d_f);
		// printVector(x, 1, N-1);

		//comp(x) = 0.5*sign(x) + 1
		evaluator->mod_switch_to_inplace(one, x.parms_id());
		evaluator->add_plain_inplace(x, one);
		evaluator->mod_switch_to_inplace(half, x.parms_id());
		evaluator->multiply_plain_inplace(x, half);
		evaluator->rescale_to_next_inplace(x);
		// cout << "multiply\n"; 
		// printVector(x, 1, N-1);
		return x;
	}

	Ciphertext initGuess(Ciphertext x)
	{
		Plaintext A, B;
		// a = 1e-3; b = 750
		// encoder->encode(-0.00019703, scale, A);
		// encoder->encode(0.14777278, scale, B);

		// a = 1e-4; b = 1000
		encoder->encode(-1.29054537e-04, scale, A);
		encoder->encode(1.29054537e-01, scale, B);

		return evalLine(x, A, B);
	}

	Ciphertext pivotTangentInitGuess(Ciphertext x)
	{
		// cout << "depth x = " << context->get_context_data(x.parms_id())->chain_index() << "\n";
		Ciphertext fact = comp(x, 7, 2);
		// cout << "depth x = " << context->get_context_data(x.parms_id())->chain_index() << "\n";
		Ciphertext l1 = evalLine(x, M1, C1);
		// cout << "depth x(l1) = " << context->get_context_data(x.parms_id())->chain_index() << "\n";
		Ciphertext l2 = evalLine(x, M2, C2);
		// cout << "depth x(l2) = " << context->get_context_data(x.parms_id())->chain_index() << "\n";

		Ciphertext res;
		Plaintext one, err2;
		encoder->encode(1.0, scale, one);
		encoder->encode(8.5e-7, scale, err);
		// encoder->encode(1.2e-6, scale, err);

		//(fact-err)*l2
		fact.scale() = scale; l1.scale() = scale; l2.scale() = scale;
		evaluator->mod_switch_to_inplace(err, fact.parms_id());
		evaluator->sub_plain_inplace(fact, err);
		// printVector(fact, 1, 5454);
		evaluator->mod_switch_to_inplace(fact, l2.parms_id());
		evaluator->multiply_inplace(l2, fact);
		evaluator->relinearize_inplace(l2, *relin_keys);
		evaluator->rescale_to_next_inplace(l2);
		// cout << "l2\n";
		
		//(1+err-fact)*l1
		evaluator->mod_switch_to_inplace(one, fact.parms_id());
		evaluator->negate_inplace(fact);
		evaluator->mod_switch_to_inplace(err, fact.parms_id());
		evaluator->add_plain_inplace(fact, err);
		evaluator->add_plain_inplace(fact, one);
		// printVector(fact, 1, 5454);
		evaluator->multiply_inplace(l1, fact);
		evaluator->relinearize_inplace(l1, *relin_keys);
		evaluator->rescale_to_next_inplace(l1);
		evaluator->add(l1, l2, res);
		// cout << "l1 " << depth(res) << "\n";
		// printVector(res, 1, 5454);
		reEncrypt(res);
		return res;
	}

	Ciphertext newtonIter(Ciphertext x, Ciphertext res, int iter=4)
	{

		for(int i=0; i < iter; i++){
			if(depth(res) < 4)
				reEncrypt(res);
			// cout << i << " " << depth(res) << "\n";
			Plaintext three_half, neg_half;
			encoder->encode(1.5, scale, three_half);
			encoder->encode(-0.5, scale, neg_half);
			
			//x^2
			Ciphertext res_sq;
			evaluator->square(res, res_sq);
			evaluator->relinearize_inplace(res_sq, *relin_keys);
			evaluator->rescale_to_next_inplace(res_sq);
			// printVector(res_sq);
			// evaluator->negate_inplace(res_sq);
			// printVector(res_sq, 3);
			// cout << "square\n";

			//-0.5*x*b
			Ciphertext res_x;
			evaluator->mod_switch_to_inplace(neg_half, x.parms_id());
			evaluator->multiply_plain(x, neg_half, res_x);
			evaluator->rescale_to_next_inplace(res_x);
			if(depth(res) < depth(res_x))
				evaluator->mod_switch_to_inplace(res_x, res.parms_id());
			else
				evaluator->mod_switch_to_inplace(res, res_x.parms_id());

			evaluator->multiply_inplace(res_x, res);
			evaluator->relinearize_inplace(res_x, *relin_keys);
			evaluator->rescale_to_next_inplace(res_x);
			// cout << "negate\n";

			//-0.5*b*x^3
			evaluator->mod_switch_to_inplace(res_sq, res_x.parms_id());
			evaluator->multiply_inplace(res_x, res_sq);
			evaluator->relinearize_inplace(res_x, *relin_keys);
			evaluator->rescale_to_next_inplace(res_x);
			// cout << "res_x\n";
			// printVector(res_x, 3);
			//1.5*x
			evaluator->mod_switch_to_inplace(three_half, res.parms_id());
			evaluator->multiply_plain_inplace(res, three_half);
			evaluator->rescale_to_next_inplace(res);			
			// cout << "constant\n";

			//-0.5*b*x^3 + 1.5*x
			evaluator->mod_switch_to_inplace(res, res_x.parms_id());
			res_x.scale() = scale; res.scale() = scale;
			evaluator->add_inplace(res, res_x);
			// cout << "final\n";
		}
		reEncrypt(res);
		return res;
	}
	pair<Ciphertext, Ciphertext> goldSchmidtIter(Ciphertext v, Ciphertext y, int d=1)
	{
		Ciphertext x, h, r, temp;
		Plaintext constant;
		encoder->encode(0.5, scale, constant);

		//GoldSchmidt's algorithm
		evaluator->mod_switch_to_inplace(y, v.parms_id());
		evaluator->multiply(v, y, x);
		evaluator->relinearize_inplace(x, *relin_keys);
		evaluator->rescale_to_next_inplace(x);
		evaluator->mod_switch_to_inplace(constant, y.parms_id());
		evaluator->multiply_plain(y, constant, h);
		evaluator->rescale_to_next_inplace(h);
		// cout << "gold\n";

		for(int i=0; i < d; i++){
			encoder->encode(0.5, scale, constant);
			//r = 0.5 - xh
			if(depth(x) < 3){
				reEncrypt(x); reEncrypt(h);
			}
			evaluator->multiply(x, h, r);
			evaluator->relinearize_inplace(r, *relin_keys);
			evaluator->rescale_to_next_inplace(r);
			r.scale() = scale;
			evaluator->negate(r, temp);
			evaluator->mod_switch_to_inplace(constant, temp.parms_id());
			evaluator->add_plain(temp, constant, r);
			// cout << "r\n";

			//x = x + x*r
			evaluator->mod_switch_to_inplace(x, r.parms_id());
			evaluator->multiply(x, r, temp);
			evaluator->relinearize_inplace(temp, *relin_keys);
			evaluator->rescale_to_next_inplace(temp);
			x.scale() = scale; temp.scale() = scale;
			evaluator->mod_switch_to_inplace(x, temp.parms_id());
			evaluator->add_inplace(x, temp);
			// cout << "x\n";

			//h = h + h*r
			evaluator->mod_switch_to_inplace(h, r.parms_id());
			evaluator->multiply(h, r, temp);
			evaluator->relinearize_inplace(temp, *relin_keys);
			evaluator->rescale_to_next_inplace(temp);
			h.scale() = scale; temp.scale() = scale;
			evaluator->mod_switch_to_inplace(h, temp.parms_id());
			evaluator->add_inplace(h, temp);			
			// cout << "h\n";
		}
		encoder->encode(2.0, scale, constant);
		evaluator->mod_switch_to_inplace(constant, h.parms_id());
		evaluator->multiply_plain_inplace(h, constant);
		evaluator->rescale_to_next_inplace(h);

		return make_pair(x, h);
	}

	Ciphertext pivotTangent(Ciphertext x, int d=8)
	{

		Ciphertext res = pivotTangentInitGuess(x);
		return newtonIter(x, res, d);
	}

	Ciphertext linearReg(Ciphertext x, int d_newt=20, int d_gold=1)
	{
		Ciphertext res = initGuess(x);
		Ciphertext y = newtonIter(x, res, d_newt);
		pair<Ciphertext, Ciphertext> sqrt_inv_sqrt = goldSchmidtIter(x, y, d_gold);
		// printVector(sqrt_inv_sqrt.first, 1);
		// printVector(sqrt_inv_sqrt.second, 1);
		return sqrt_inv_sqrt.second;
	}

};

int main()
{
	//Initialisation parameters
	EncryptionParameters parms(scheme_type::ckks);
	long logN = 14;
	size_t poly_modulus_degree = 1 << logN;
	double scale = pow(2.0, 40);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 58, 40, 40, 40, 40, 40, 40, 40, 40, 58 }));
    SEALContext context(parms);

	KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    CKKSEncoder encoder(context);		
    RelinKeys relin_keys;
	keygen.create_relin_keys(relin_keys);
	// GaloisKeys galois_keys;
    // keygen.create_galois_keys(galois_keys);

    InvSqrt inv_sqrt(context, encryptor, decryptor, encoder, evaluator, scale, relin_keys);

    size_t slot_count = encoder.slot_count();// = poly_modulus_degree/2
    double a = 1e-4, b=1e3;
    // double a = 1e-3, b=750;
    vector<double> c(slot_count);
	double h1 = (1-a)/(slot_count/2-1), h2 = (b-1)/(slot_count/2-1);
    for(int i=0; i < slot_count/2; i++){
    	c[i] = a + i*h1;
    	c[i+slot_count/2] = 1 + i*h2;
    	// cout << c[i] << " " << c[i+slot_count/2] << "\n";
    }

    Plaintext plain;
    encoder.encode(c, scale, plain);
    Ciphertext x, res;
    encryptor.encrypt(plain, x);

    vector<double> y(slot_count), y1(slot_count), y2(slot_count), l(slot_count);
    res = inv_sqrt.linearReg(x, 15, 5);
    decryptor.decrypt(res, plain);
    encoder.decode(plain, y2);

    res = inv_sqrt.pivotTangent(x, 9);
    decryptor.decrypt(res, plain);
    encoder.decode(plain, y1);

	res = inv_sqrt.pivotTangentInitGuess(x);
    decryptor.decrypt(res, plain);
    encoder.decode(plain, l);

    double err1 = 0, err2 = 0;
    for(int i=0; i < slot_count; i++){
    	// cout << c[i] << " " << y1[i] << " " << y2[i] << " " << pow(c[i], -0.5) << "\n";
    	y[i] = pow(c[i], -0.5);
    	err1 += abs(pow(c[i], -0.5)-y1[i]);
    	err2 += abs(pow(c[i], -0.5)-y2[i]);
    }
    double k1 = inv_sqrt.getBounds(1), k2 = inv_sqrt.getBounds(2);

    cout << "Pivot-Tangent Method error = " << err1/slot_count << "\n";
    cout << "Linear Regression Method error = " << err2/slot_count << "\n";
    // cout << y1[0] << " " << y1[slot_count-1] << " " << y2[0] << " " << y2[slot_count-1] << "\n";

 	plt::named_plot("Original", c, y, "go");
 	plt::named_plot("Predicted", c, y1);
 	plt::named_plot("Init Guess", c, l);
 	plt::plot(c, [&](double d){return k1*pow(d, -0.5);}, "r--");
 	plt::plot(c, [&](double d){return k2*pow(d, -0.5);}, "r--");
 	plt::xlim(-0.5, 8.0);
 	plt::ylim(-0.5, 8.0);
 	plt::legend();
 	plt::show();   
}