// lib includes
#include "wxmpimsg.h"

WxMpiMsg::WxMpiMsg()
  : WxMsgBase(0, MPI_ANY_TAG), _comm(MPI_COMM_WORLD) {
  // add templated messengers to the base class
  this->addMsg<char>
    ( new WxMpiMsgTmpl<char>(_comm) );
  this->addMsg<unsigned char>
    ( new WxMpiMsgTmpl<unsigned char>(_comm) );
  this->addMsg<short>
    ( new WxMpiMsgTmpl<short>(_comm) );
  this->addMsg<unsigned short>
    ( new WxMpiMsgTmpl<unsigned short>(_comm) );
  this->addMsg<int>
    ( new WxMpiMsgTmpl<int>(_comm) );
  this->addMsg<unsigned int>
    ( new WxMpiMsgTmpl<unsigned int>(_comm) );
  this->addMsg<long>
    ( new WxMpiMsgTmpl<long>(_comm) );
  this->addMsg<unsigned long>
    ( new WxMpiMsgTmpl<unsigned long>(_comm) );
  //this->addMsg<float>
  //  ( new WxMpiMsgTmpl<float>(_comm) );
  this->addMsg<double>
    ( new WxMpiMsgTmpl<double>(_comm) );
  this->addMsg<long double>
    ( new WxMpiMsgTmpl<long double>(_comm) );
  this->addMsg<long long int>
    ( new WxMpiMsgTmpl<long long int>(_comm) );
}

WxMpiMsg::WxMpiMsg(WxMpiMsg *parent, MPI_Comm comm)
  : WxMsgBase(0, MPI_ANY_TAG, parent) {
  // add templated messengers to the base class
  this->addMsg<char>
    ( new WxMpiMsgTmpl<char>(_comm) );
  this->addMsg<unsigned char>
    ( new WxMpiMsgTmpl<unsigned char>(_comm) );
  this->addMsg<short>
    ( new WxMpiMsgTmpl<short>(_comm) );
  this->addMsg<unsigned short>
    ( new WxMpiMsgTmpl<unsigned short>(_comm) );
  this->addMsg<int>
    ( new WxMpiMsgTmpl<int>(_comm) );
  this->addMsg<unsigned int>
    ( new WxMpiMsgTmpl<unsigned int>(_comm) );
  this->addMsg<long>
    ( new WxMpiMsgTmpl<long>(_comm) );
  this->addMsg<unsigned long>
    ( new WxMpiMsgTmpl<unsigned long>(_comm) );
  //this->addMsg<float>
  //  ( new WxMpiMsgTmpl<float>(_comm) );
  this->addMsg<double>
    ( new WxMpiMsgTmpl<double>(_comm) );
  this->addMsg<long double>
    ( new WxMpiMsgTmpl<long double>(_comm) );
  this->addMsg<long long int>
    ( new WxMpiMsgTmpl<long long int>(_comm) );
}

WxMsgBase*
WxMpiMsg::createSubComm(const std::vector<int>& ranks) 
{
  MPI_Group old_group, new_group;
  MPI_Comm new_comm;
  int *my_ranks = new int[ranks.size()];
  for (unsigned i=0; i<ranks.size(); ++i)
    my_ranks[i] = ranks[i];
  MPI_Comm_group(_comm, &old_group);
  MPI_Group_incl(old_group, ranks.size(), my_ranks, &new_group);
  MPI_Comm_create(_comm, new_group, &new_comm);
  delete [] my_ranks;

  if (new_comm == MPI_COMM_NULL)
    return 0;
  WxMsgBase *c = new WxMpiMsg(this, new_comm);
  return c;
}

void * 
WxMpiMsg::finishRecv(WxMsgStatus wxms)
{
  WxMpiMsgStatus_v *ms = static_cast<WxMpiMsgStatus_v *>(wxms);
  MPI_Status status;
  MPI_Wait(&ms->request, &status);
  void *data = ms->data;
  delete wxms;
  return data;
}

bool
WxMpiMsg::checkRecv(WxMsgStatus wxms)
{
  WxMpiMsgStatus_v *ms = static_cast<WxMpiMsgStatus_v *>(wxms);
  MPI_Status status;
  int flag;
  MPI_Test(&ms->request, &flag, &status);
  return (flag != 0) ? true : false;
}
