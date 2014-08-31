#ifndef __wxsubsolver__h__
#define __wxsubsolver__h__

// WarpX lib includes
#include <wxcrypt.h>
#include <wxmsgbase.h>
#include <wxobject.h>
#include <wxstepper.h>

// std includes
#include <string>
#include <typeinfo>
#include <vector>

// forward declare solver
template <typename REAL> class ApSolver;

/**
 * Base class for sub-solvers in WarpX system. A subsolver is an
 * abject which encodes algorithms to transform input variables into
 * output variables. It is like a subroutine in a programming
 * langauge.
 */
template <typename REAL>
class WxSubSolver : public WxStepper<REAL>
{
  public:
/**
 * Create new subsolver with given name
 *
 * @param name Name of subsolver
 */
    WxSubSolver(const std::string& name);

/** Destroy subsolver */
    virtual ~WxSubSolver();

/** Return type of subsolver */
    virtual const std::type_info& type() const;

/**
 * Setup subsolver object using supplied cryptset
 *
 * @param wxc Cryptset to use for setting
 */
    virtual void setup(const WxCryptSet& wxc);

/**
 * Declare the types of read and write variables accepted by this
 * subsolver
 */
    virtual void declareTypes();

/**
 * Get names of all read variables used
 *
 * @return names of read variables
 */
    std::vector<std::string> readVarNames() const;

/**
 * Get names of all write variables used
 *
 * @return names of write variables
 */
    std::vector<std::string> writeVarNames() const;

/**
 * Set type information for a read variable
 *
 * @param pos Position of read variable
 * @param type its type gotten from typeid operator
 */
    void setReadVarType(unsigned pos, const std::type_info& type);

/**
 * Set type information for a read variable
 *
 * @param pos Position of read variable
 * @return type of read variable
 */
    const std::type_info& getReadVarType(unsigned pos);

/**
 * Set type information for a write variable
 *
 * @param pos Position of write variable
 * @param type its type gotten from typeid operator
 */
    void setWriteVarType(unsigned pos, const std::type_info& type);

/**
 * Set type information for a write variable
 *
 * @param pos Position of write variable
 * @return type of write variable
 */
    const std::type_info& getWriteVarType(unsigned pos);

/**
 * If the subsolver takes an arbitrary number of read variables then
 * set type information for the variables following the last required
 * read variable.
 *
 * @param type its type gotten from typeid operator
 */
    void setLastReadVarType(const std::type_info& type);

/**
 * If the subsolver takes an arbitrary number of read variables then
 * return type of last positional variables
 *
 * @return type of last read variables
 */
    const std::type_info& getLastReadVarType();

/**
 * If the subsolver takes an arbitrary number of write variables then
 * set type information for the variables following the last required
 * write variable.
 *
 * @param type its type gotten from typeid operator
 */
    void setLastWriteVarType(const std::type_info& type);

/**
 * If the subsolver takes an arbitrary number of write variables then
 * return type of last positional variables
 *
 * @return type of last write variables
 */
    const std::type_info& getLastWriteVarType();

/**
 * Return number of read variables expected by this subsolver
 *
 * @return number of read variables expected
 */
    unsigned numReadVars() const;

/**
 * Return number of write variables expected by this subsolver
 *
 * @return number of write variables expected
 */
    unsigned numWriteVars() const;

/**
 * Does this subsolver have variable number of read variables?
 *
 * @return true if variable number of read variables
 */
    bool hasVariableReadVars() const;

/**
 * Does this subsolver have variable number of write variables?
 *
 * @return true if variable number of write variables
 */
    bool hasVariableWriteVars() const;

/**
 * Return true if this updater needs a grid to work on
 *
 * @return true if a grid is needed, false otherwise.
 */
    virtual bool needsGrid() const;

/**
 * Set parent solver object
 * 
 * @param parent pointer to parent
 */
    void setParent(ApSolver<REAL> *parent);

/**
 * Get parent solver object
 * 
 * @return pointer to parent
 */
    ApSolver<REAL>* getParent() const;

/**
 * Get number of actual read variables used when user calls this
 * subsolver
 *
 * @param number of read variables
 */
    unsigned numActualReadVars() const;

/**
 * Get number of actual write variables used when user calls this
 * subsolver
 *
 * @param number of write variables
 */
    unsigned numActualWriteVars() const;

/**
 * Get grid on which solver works
 *
 * @return domain name
 */
    //WxGridBox<REAL> getGrid() const;

/**
 * Returns list of write variables for this subsolver
 *
 * @return List of write variables
 */
    std::vector<std::string> getWriteVarList() const;

/**
 * Returns list of read variables for this subsolver
 *
 * @return List of read variables
 */
    std::vector<std::string> getReadVarList() const;

/**
 * Sets list of write variables for this subsolver
 *
 * @param wr List of write vairables
 */
    void setWriteVarList(const std::vector<std::string>& wr);

/**
 * Sets list of read variables for this subsolver
 *
 * @param wr List of write vairables
 */
    void setReadVarList(const std::vector<std::string>& rv);

  protected:

/**
 * Get a reference to a read variable
 *
 * @param pos position of variable in readVars list
 * @return reference to variable
 */
    template <typename UT>
    const UT& getReadVar(unsigned pos) const {
      return this->_parent->template getConstVar<UT>(this->_readVars[pos]);
    }

/**
 * Get a reference to a write variable
 *
 * @param pos position of variable in writeVars list
 * @return reference to variable
 */
    template <typename UT>
    UT& getWriteVar(unsigned pos) const {
      return this->_parent->template getVar<UT>(this->_writeVars[pos]);
    }

  private:
    ApSolver<REAL> *_parent; // parent solver
    std::vector<std::string> _readVars; // variables to read from
    std::vector<std::string> _writeVars; // variables to write to
    std::map<unsigned, const std::type_info*> _readVarType; // map of read variables to type
    std::map<unsigned, const std::type_info*> _writeVarType; // map of write variables to type
    const std::type_info* _lastReadVarType; // type of last read variable
    const std::type_info* _lastWriteVarType; // type of last write variable
    std::string _onGrid; // grid on which solver acts
};

#endif //  __wxsubsolver__h__
