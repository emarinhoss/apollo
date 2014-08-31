#ifndef __wxmsgop__
#define __wxmsgop__

// list of supported operations for all-reduce operation
enum WxMsgOp
{
    WX_MSG_NOP,
    WX_MSG_MIN,
    WX_MSG_MAX,
    WX_MSG_AND,
    WX_MSG_SUM
};

/**
 * Generic reduction operation wrapper
 */
class WxMsgOpC
{
  protected:
/**
 * Constructor
 */
    WxMsgOpC() {}
};

#endif //  __wxmsgop__
