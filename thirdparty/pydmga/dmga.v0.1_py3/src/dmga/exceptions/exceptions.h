/* 
 * File:   exceptions.h
 * Author: Robson
 * 
 * This file contains all exception classes that can be thrown from the routines 
 * in DMGa library
 *
 * Created on 16 kwiecień 2012, 11:30
 */

#ifndef EXCEPTIONS_H
#define	EXCEPTIONS_H

#include <ostream>

namespace dmga{

namespace exceptions{

class BaseException{	
public:
	std::string what;
	BaseException* cause;

	BaseException(const std::string& what, BaseException& cause) : what(what), cause(cause.makeCause()) { //this is not nice (this pointers), but should work due to the fact how the stack untangling after exception is done
	}

	BaseException(const std::string& what = "") : what(what), cause(0){
	}

	virtual BaseException* makeCause(){
		if (cause){
			return new BaseException(what, *cause);
		}else{
			return new BaseException(what);
		}
	}

	virtual ~BaseException(){
		if (cause){
			delete cause;
		}
	}

	friend std::ostream& operator<<(std::ostream& out, BaseException& exception){
		out << exception.what << std::endl;
		if (exception.cause){
			out << "\tCaused by " << (*exception.cause);
		}
		return out;
	}
};

class NotImplementedYet : public BaseException{
public:
	NotImplementedYet(const std::string& what = "") : BaseException("NotImplementedYetException: " + what)
	{}
	NotImplementedYet(const std::string& what, BaseException& cause) : BaseException("NotImplementedYetException: " + what, cause)
	{}
};

class InvalidInput : public BaseException{
public:
	InvalidInput(const std::string& what = "") : BaseException("InvalidInputException: " + what)
	{}
	InvalidInput(const std::string& what, BaseException& cause) : BaseException("InvalidInputException: " + what, cause)
	{}
};

class NotConfigured : public BaseException{
public:
	NotConfigured(const std::string& what = "") : BaseException("NotConfiguredException: " + what)
	{}
	NotConfigured(const std::string& what, BaseException& cause) : BaseException("NotConfiguredException: " + what, cause)
	{}
};

class LogicViolation : public BaseException{
public:
	LogicViolation(const std::string& what = "") : BaseException("ApplicationLogicException: " + what)
	{}
	LogicViolation(const std::string& what, BaseException& cause) : BaseException("ApplicationLogicException: " + what, cause)
	{}
};

class NotSupported : public BaseException{
public:
	NotSupported(const std::string& what = "") : BaseException("NotSupportedException: " + what)
	{}
	NotSupported(const std::string& what, BaseException& cause) : BaseException("NotSupportedException: " + what, cause)
	{}
};

/**
 * this is for similar concept as in python iterators
 */
class StopIteration : public BaseException{
public:
	StopIteration(const std::string& what = "") : BaseException("StopIteration: " + what)
	{}
	StopIteration(const std::string& what, BaseException& cause) : BaseException("NotSupportedException: " + what, cause)
	{}
};

} //namespace exceptions

} //namespace dmga


#endif	/* EXCEPTIONS_H */

