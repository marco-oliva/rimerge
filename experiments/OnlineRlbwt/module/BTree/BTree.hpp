/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file BTree.hpp
 * @brief Pointer-based implementation of upper part of B+tree.
 * @attention Bottom part of B+tree can be implemented in more space efficient way (e.g., array-based implementation).
 * @author Tomohiro I
 * @date 2017-01-26
 * @todo Support delete.
 */
#ifndef INCLUDE_GUARD_BTree
#define INCLUDE_GUARD_BTree

#include <algorithm>
#include <cassert>
#include <iostream>

namespace itmmti
{
  template<class BTreeNodeT>
  struct SuperRoot
  {
    BTreeNodeT * root_;


    SuperRoot() : root_(nullptr)
    {}


    void setRoot
    (
     BTreeNodeT * root
     ) noexcept {
      assert(root->isRoot());
      assert(root->isUnderSuperRoot());

      root_ = root;
      root->setParentRef(reinterpret_cast<BTreeNodeT *>(this), root->getIdxInSibling());
    }
  };



  
  /*!
   * @brief Pointer-based implementation of upper part of B+tree.
   * @attention Bottom part of B+tree can be implemented in more space efficient way (e.g., using packed array to store weights).
   * @tparam kB should be in {4, 8, 16, 32, 64, 128}. kB/2 <= 'numChildren_' <= kB.
   */
  template <uint8_t tparam_kB = 64>
  class BTreeNode
  {
  public:
    //// Public constant, alias etc.
    static constexpr uint8_t kB{tparam_kB};
    using BTreeNodeT = BTreeNode<kB>;
    using SuperRootT = SuperRoot<BTreeNodeT>;
    static constexpr uintptr_t NOTFOUND{UINTPTR_MAX};


  private:
    //// Private constant, alias etc.
    /*!
     * @brief For representing current status of node by 'flags_'.
     */
    enum {
      isRootBit = 1, //!< This node is a root.
      isBorderBit = 2, //!< This node is on border, i.e., its children are bottom nodes.
      isJumpToBtmBit = 4, //!< If true, "lmJumpNode_" points to leftmost bottom node. Otherwise, it points to leftmost border node.
      isUnderSuperRootBit = 8, //!< This node is in a subtree whose root is a child of super root.
      isDummyBit = 16, //!< This node is a dummy node.
    };


  private:
    //// Private member variables
    uint64_t psum_[kB+1]; //!< Partial sum: psum_[i+1] = sum_{i = 0}^{i} [weight of i-th child (0base)]
    BTreeNodeT * parent_; //!< Pointer to parent node.
    uint8_t idxInSibling_; //!< This node is 'idxInSibling_'-th child (0base) of its parent.
    uint8_t numChildren_; //!< Current num of children.
    uint8_t flags_;
    BTreeNodeT * children_[kB]; //!< Pointers to children. Note that the children might not BTreeNode type when this node is in border.
    BTreeNodeT * lmJumpNode_; //!< To remember leftmost border node of this node.


  public:
    //// Constructor
    BTreeNode
    (
     void * lmJumpNode,
     bool isRoot,
     bool isBorder,
     bool isJumpToBtm,
     bool isUnderSuperRoot,
     bool isDummy = false
     )
    : parent_(nullptr),
      idxInSibling_(0),
      numChildren_(0),
      flags_(isRoot * isRootBit | isBorder * isBorderBit | isJumpToBtm * isJumpToBtmBit | isUnderSuperRoot * isUnderSuperRootBit | isDummy * isDummyBit),
      lmJumpNode_(reinterpret_cast<BTreeNodeT *>(lmJumpNode))
    {
      psum_[0] = 0;
      if (isBorder && !isJumpToBtm) {
        lmJumpNode_ = this;
      }
    }
    ~BTreeNode() = default;
    BTreeNode(const BTreeNode &) = delete;
    BTreeNode & operator=(const BTreeNode &) = delete;


    /*!
     * @brief Clear all upper nodes in B+tree.
     * @node Bottom nodes might be cleared manually before calling this function.
     */
    void clearUpperNodes() noexcept {
      if (!isBorder()) {
        for (uint8_t i = 0; i < numChildren_; ++i) {
          children_[i]->clearUpperNodes();
        }
      }
      delete this;
    }


  public:
    //// simple getter
    /*!
     * @brief Get the sum of weights of child subtrees numbered from 0 to i-1 (note that i is NOT included).
     */
    uint64_t getPSum
    (
     const uint8_t idx_excl //!< in [0.."numChildren_"].
     ) const noexcept {
      assert(idx_excl <= numChildren_);

      return psum_[idx_excl];
    }


    /*!
     * @brief Get the weight of a child subtree.
     */
    uint64_t getWeightOfChild
    (
     const uint8_t idx //!< in [0.."numChildren_").
     ) const noexcept {
      assert(idx < numChildren_);

      return psum_[idx+1] - psum_[idx];
    }


    /*!
     * @brief Get const psum array.
     */
    const uint64_t * getConstPtr_psum () const noexcept {
      return psum_;
    }


    /*!
     * @brief Get the weight of this node.
     */
    uint64_t getSumOfWeight() const noexcept {
      return psum_[numChildren_];
    }


    /*!
     * @brief Get pointer to the i-th child (0base).
     */
    BTreeNodeT * getChildPtr
    (
     uint8_t i //!< in [0..numChildren_).
     ) const noexcept {
      assert(i < numChildren_);

      return children_[i];
    }


    /*!
     * @brief Get pointer to parent.
     */
    BTreeNodeT * getParent() const noexcept {
      return parent_;
    }


    /*!
     * @brief Get idxInSibling_.
     */
    uint8_t getIdxInSibling() const noexcept {
      return idxInSibling_;
    }


    /*!
     * @brief Get num of children.
     */
    uint8_t getNumChildren() const noexcept {
      return numChildren_;
    }


    /*!
     * @brief Return if this node is root.
     */
    bool isRoot() const noexcept {
      return flags_ & isRootBit;
    }


    /*!
     * @brief Return if this node is in border.
     */
    bool isBorder() const noexcept {
      return flags_ & isBorderBit;
    }


    /*!
     * @brief Return if "lmJumpNode_" of this node points to leftmost bottom node.
     */
    bool isJumpToBtm() const noexcept {
      return flags_ & isJumpToBtmBit;
    }


    /*!
     * @brief Return if this node is in subtree whose root is a child of super root.
     */
    bool isUnderSuperRoot() const noexcept {
      return flags_ & isUnderSuperRootBit;
    }


    /*!
     * @brief Return if this node is a dummy node.
     */
    bool isDummy() const noexcept {
      return flags_ & isDummyBit;
    }


  public:
    //////////////////////////////// functions to traverse tree
    /*!
     * @brief Get leftmost jump node.
     * @node If "isJumpToBtm() == true", it points to leftmost bottom node. Otherwise, it points to leftmost border node.
     */
    const BTreeNodeT * getLmJumpNode() const noexcept {
      return lmJumpNode_;
    }


    BTreeNodeT * getLmJumpNode() noexcept {
      return lmJumpNode_;
    }


    /*!
     * @brief Get leftmost border node. It takes O(1) time.
     */
    const BTreeNodeT * getLmBorderNode_DirectJump() const noexcept {
      assert(!isJumpToBtm());

      return lmJumpNode_;
    }


    BTreeNodeT * getLmBorderNode_DirectJump() noexcept {
      assert(!isJumpToBtm());

      return lmJumpNode_;
    }


    /*!
     * @brief Get leftmost border node. It takes O(h) time, where h is the height of this node.
     */
    const BTreeNodeT * getLmBorderNode_NoDirectJump() const noexcept {
      const auto * node = this;
      while (!(node->isBorder())) {
        node = node->children_[0];
      }
      return node;
    }


    BTreeNodeT * getLmBorderNode_NoDirectJump() noexcept {
      return const_cast<BTreeNodeT *>(static_cast<const BTreeNodeT &>(*this).getLmBorderNode_NoDirectJump());
    }


    /*!
     * @brief Get rightmost border node. It takes O(h) time, where h is the height of this node.
     */
    const BTreeNodeT * getRmBorderNode() const noexcept {
      const auto * node = this;
      while (!(node->isBorder())) {
        node = node->children_[node->numChildren_ - 1];
      }
      return node;
    }


    BTreeNodeT * getRmBorderNode() noexcept {
      return const_cast<BTreeNodeT *>(static_cast<const BTreeNodeT &>(*this).getRmBorderNode());
    }


    /*!
     * @brief Get leftmost bottom. It takes O(1) time.
     */
    const BTreeNodeT * getLmBtm_DirectJump() const noexcept {
      assert(isJumpToBtm());

      return lmJumpNode_;
    }


    BTreeNodeT * getLmBtm_DirectJump() noexcept {
      assert(isJumpToBtm());

      return lmJumpNode_;
    }


    /*!
     * @brief Get leftmost bottom. It takes O(1) time.
     */
    const BTreeNodeT * getLmBtm_NoDirectJump() const noexcept {
      assert(!isJumpToBtm());

      return lmJumpNode_->getChildPtr(0);
    }


    BTreeNodeT * getLmBtm_NoDirectJump() noexcept {
      assert(!isJumpToBtm());

      return lmJumpNode_->getChildPtr(0);
    }


    /*!
     * @brief Get rightmost bottom. It takes O(h) time, where h is the height of this node.
     */
    const BTreeNodeT * getRmBtm() const noexcept {
      const auto * node = this;
      while (!(node->isBorder())) {
        node = node->children_[node->getNumChildren() - 1];
      }
      return node->children_[node->getNumChildren() - 1];
    }


    BTreeNodeT * getRmBtm() noexcept {
      return const_cast<BTreeNodeT *>(static_cast<const BTreeNodeT &>(*this).getRmBtm());
    }


    /*!
     * @brief Return next bottm starting from "idxInSib"-th child (0base) of this node.
     */
    const BTreeNodeT * getNextBtm_DirectJump
    (
     uint8_t idxInSib
     ) const noexcept {
      assert(isJumpToBtm());

      const auto * node = this;
      while (idxInSib + 1 == node->getNumChildren() && !(node->isRoot())) {
        idxInSib = node->getIdxInSibling();
        node = node->getParent();
      }
      if (idxInSib + 1 < node->getNumChildren()) {
        if (node->isBorder()) {
          return node->getChildPtr(idxInSib + 1);
        } else {
          return node->getChildPtr(idxInSib + 1)->getLmBtm_DirectJump();
        }
      }
      return reinterpret_cast<BTreeNodeT *>(NOTFOUND);
    }


    BTreeNodeT * getNextBtm_DirectJump
    (
     uint8_t idxInSib
     ) noexcept {
      return const_cast<BTreeNodeT *>(static_cast<const BTreeNodeT &>(*this).getNextBtm_DirectJump(idxInSib));
    }


    /*!
     * @brief Return next bottm starting from "idxInSib"-th child (0base) of this node.
     */
    const BTreeNodeT * getNextBtm_NoDirectJump
    (
     uint8_t idxInSib
     ) const noexcept {
      assert(!isJumpToBtm());

      const auto * node = this;
      while (idxInSib + 1 == node->getNumChildren() && !(node->isRoot())) {
        idxInSib = node->getIdxInSibling();
        node = node->getParent();
      }
      if (idxInSib + 1 < node->getNumChildren()) {
        if (node->isBorder()) {
          return node->getChildPtr(idxInSib + 1);
        } else {
          return node->getChildPtr(idxInSib + 1)->getLmBtm_NoDirectJump();
        }
      }
      return reinterpret_cast<BTreeNodeT *>(NOTFOUND);
    }


    BTreeNodeT * getNextBtm_NoDirectJump
    (
     uint8_t idxInSib
     ) noexcept {
      return const_cast<BTreeNodeT *>(static_cast<const BTreeNodeT &>(*this).getNextBtm_NoDirectJump(idxInSib));
    }


    /*!
     * @brief Return previous bottm starting from "idxInSib"-th child (0base) of this node.
     */
    const BTreeNodeT * getPrevBtm
    (
     uint8_t idxInSib
     ) const noexcept {
      const auto * node = this;
      while (idxInSib == 0 && !(node->isRoot())) {
        idxInSib = node->getIdxInSibling();
        node = node->getParent();
      }
      if (idxInSib) {
        if (node->isBorder()) {
          return node->getChildPtr(idxInSib - 1);
        } else {
          return node->getChildPtr(idxInSib - 1)->getRmBtm();
        }
      }
      return reinterpret_cast<BTreeNodeT *>(NOTFOUND);
    }


    BTreeNodeT * getPrevBtm
    (
     uint8_t idxInSib
     ) noexcept {
      return const_cast<BTreeNodeT *>(static_cast<const BTreeNodeT &>(*this).getPrevBtm(idxInSib));
    }


    /*!
     * @brief Return next bottm starting from "idxInSib"-th child (0base) of this node.
     */
    const BTreeNodeT * getNextBtmRef_DirectJump
    (
     uint8_t & idxInSib //!< [in,out]
     ) const noexcept {
      assert(!isJumpToBtm());

      auto * node = this;
      while (++idxInSib == node->getNumChildren() && !(node->isRoot())) {
        idxInSib = node->getIdxInSibling();
        node = node->getParent();
      }
      if (idxInSib < node->getNumChildren()) {
        if (!(node->isBorder())) {
          node = node->getChildPtr(idxInSib);
          node = node->getLmBorderNode_DirectJump();
          idxInSib = 0;
        }
        return node;
      }
      return reinterpret_cast<BTreeNodeT *>(NOTFOUND);
    }


    BTreeNodeT * getNextBtmRef_DirectJump
    (
     uint8_t & idxInSib //!< [in,out]
     ) noexcept {
      return const_cast<BTreeNodeT *>(static_cast<const BTreeNodeT &>(*this).getNextBtmRef_DirectJump(idxInSib));
    }


    /*!
     * @brief Return next bottm starting from "idxInSib"-th child (0base) of this node.
     */
    const BTreeNodeT * getNextBtmRef_NoDirectJump
    (
     uint8_t & idxInSib //!< [in,out]
     ) const noexcept {
      assert(isJumpToBtm());

      auto * node = this;
      while (++idxInSib == node->getNumChildren() && !(node->isRoot())) {
        idxInSib = node->getIdxInSibling();
        node = node->getParent();
      }
      if (idxInSib < node->getNumChildren()) {
        if (!(node->isBorder())) {
          node = node->getChildPtr(idxInSib);
          node = node->getLmBorderNode_NoDirectJump();
          idxInSib = 0;
        }
        return node;
      }
      return reinterpret_cast<BTreeNodeT *>(NOTFOUND);
    }


    BTreeNodeT * getNextBtmRef_NoDirectJump
    (
     uint8_t & idxInSib //!< [in,out]
     ) noexcept {
      return const_cast<BTreeNodeT *>(static_cast<const BTreeNodeT &>(*this).getNextBtmRef_NoDirectJump(idxInSib));
    }


    /*!
     * @brief Return previous bottm starting from "idxInSib"-th child (0base) of this node.
     */
    const BTreeNodeT * getPrevBtmRef
    (
     uint8_t & idxInSib //!< [in,out]
     ) const noexcept {
      auto * node = this;
      while (idxInSib == 0 && !(node->isRoot())) {
        idxInSib = node->getIdxInSibling();
        node = node->getParent();
      }
      if (idxInSib--) {
        if (!isBorder()) {
          node = node->getChildPtr(idxInSib);
          node = node->getRmBorderNode();
          idxInSib = node->getNumChildren();
        }
        return node;
      }
      return reinterpret_cast<BTreeNodeT *>(NOTFOUND);
    }


    BTreeNodeT * getPrevBtmRef
    (
     uint8_t & idxInSib //!< [in,out]
     ) noexcept {
      return const_cast<BTreeNodeT *>(static_cast<const BTreeNodeT &>(*this).getPrevBtmRef(idxInSib));
    }


    /*!
     * @brief Return next sibling of this node.
     */
    const BTreeNodeT * getNextSib() const noexcept {
      uint32_t level = 1;
      const auto * node = this;
      while (!(node->isRoot())) {
        const auto idxInSib = node->getIdxInSibling();
        node = node->getParent();
        if (idxInSib + 1 < node->getNumChildren()) {
          node = node->getChildPtr(idxInSib + 1);
          while (--level) {
            node = node->getChildPtr(0);
          }
          return node;
        }
        ++level;
      }
      return reinterpret_cast<BTreeNodeT *>(NOTFOUND);
    }


    /*!
     * @brief Return next sibling of this node.
     */
    BTreeNodeT * getNextSib() noexcept {
      return const_cast<BTreeNodeT *>(static_cast<const BTreeNodeT &>(*this).getNextSib());
    }


    /*!
     * @brief Return previous sibling of this node.
     */
    const BTreeNodeT * getPrevSib () const noexcept {
      uint32_t level = 1;
      const auto * node = this;
      while (!(node->isRoot())) {
        const auto idxInSib = node->getIdxInSibling();
        node = node->getParent();
        if (idxInSib) {
          node = node->getChildPtr(idxInSib - 1);
          while (--level) {
            node = node->getChildPtr(node->getNumChildren() - 1);
          }
          return node;
        }
        ++level;
      }
      return reinterpret_cast<BTreeNodeT *>(NOTFOUND);
    }


    /*!
     * @brief Return previous sibling of this node.
     */
    BTreeNodeT * getPrevSib () noexcept {
      return const_cast<BTreeNodeT *>(static_cast<const BTreeNodeT &>(*this).getPrevSib());
    }


    /*!
     * @brief
     *   Return partial sum (or row number "ROW") up to the node
     *   (inclusive iff "inclusive == true") indicated by "idx"-th child (0base) of this node.
     */
    /*
      template<uint8_t ROW>
      uint64_t calcPSum
      (
      uint8_t idx,
      bool inclusive = false
      ) const noexcept {
      static_assert(ROW < ROW_NUM, "ROW should be smaller than ROW_NUM");
      assert(isBorder());

      const auto * node = this;
      uint64_t ret = 0;
      while (true) {
      ret += node->getPSum<ROW>(idx + inclusive);
      if (node->isRoot()) {
      return ret;
      }
      idx = node->getIdxInSibling();
      node = node->getParent();
      }
      }
    */


    /*!
     * @brief Return partial sum up to the node (exclusive) indicated by "idx"-th child (0base) of this node.
     */
    uint64_t calcPSum
    (
     uint8_t idx
     ) const noexcept {
      assert(isBorder());

      const auto * node = this;
      uint64_t ret = 0;
      while (true) {
        ret += node->getPSum(idx);
        if (node->isRoot()) {
          return ret;
        }
        idx = node->getIdxInSibling();
        node = node->getParent();
      }
    }


    /*!
     * @brief Return partial sum up to the node (exclusive) indicated by "idx"-th child (0base) of this node.
     */
    uint64_t calcPSum
    (
     uint8_t idx,
     const BTreeNodeT *& retNode //!< [out] To capture the root of the BTree.
     ) const noexcept {
      assert(isBorder());

      retNode = this;
      uint64_t ret = 0;
      while (true) {
        ret += retNode->getPSum(idx);
        if (retNode->isRoot()) {
          return ret;
        }
        idx = retNode->getIdxInSibling();
        retNode = retNode->getParent();
      }
    }


    /*!
     * @brief Return partial sum up to the node (exclusive) indicated by "idx"-th child (0base) of this node.
     */
    uint64_t calcPSum
    (
     uint8_t idx,
     BTreeNodeT *& retNode //!< [out] To capture the root of the BTree.
     ) noexcept {
      assert(isBorder());

      retNode = this;
      uint64_t ret = 0;
      while (true) {
        ret += retNode->getPSum(idx);
        if (retNode->isRoot()) {
          return ret;
        }
        idx = retNode->getIdxInSibling();
        retNode = retNode->getParent();
      }
    }


    /*!
     * @brief Traverse tree looking for "pos" in weights array.
     * @return Pointer to bottom node where partial sum of "pos" is achieved, where weight-0 nodes (e.g. dummy nodes) are skipped.
     */
    const BTreeNodeT * searchPos
    (
     uint64_t & pos //!< [in,out] Give global position to search. It is modified to relative position in bottom node.
     ) const noexcept {
      assert(pos < this->getSumOfWeight());

      const BTreeNodeT * node = this;
      while (true) {
        uint8_t i = 0;
        auto array = node->getConstPtr_psum();
        while (pos >= array[i + 1]) {
          ++i;
        }
        pos -= array[i];
        if (node->isBorder()) {
          return node->getChildPtr(i);
        }
        node = node->getChildPtr(i);
      }
    }


    /*!
     * @brief Traverse tree looking for "pos" in weights array.
     * @return Pointer to bottom node where partial sum of "pos" is achieved, where weight-0 nodes (e.g. dummy nodes) are skipped.
     */
    BTreeNodeT * searchPos
    (
     uint64_t & pos //!< [in,out] Give global position to search. It is modified to relative position in bottom node.
     ) noexcept {
      return const_cast<BTreeNodeT *>(static_cast<const BTreeNodeT &>(*this).searchPos(pos));
    }


    /*!
     * @brief Traverse tree looking for "pos" in weights array.
     * @return Pointer to bottom node where partial sum of "pos" is achieved, where weight-0 nodes (e.g. dummy nodes) are skipped.
     */
    void searchPos
    (
     uint64_t & pos, //!< [in,out] Give global position to search. It is modified to relative position in bottom node.
     const BTreeNodeT *& retNode, //!< [out] To capture parent of returned bottom node.
     uint8_t & retIdx //!< [out] To capture sibling idx of returned bottom node.
     ) const noexcept {
      assert(pos < this->getSumOfWeight());

      retNode = this;
      while (true) {
        retIdx = 0;
        auto array = retNode->getConstPtr_psum();
        while (pos >= array[retIdx + 1]) {
          ++retIdx;
        }
        pos -= array[retIdx];
        if (retNode->isBorder()) {
          return;
        }
        retNode = retNode->getChildPtr(retIdx);
      }
    }


    /*!
     * @brief Traverse tree looking for "pos" in weights array.
     * @return Pointer to bottom node where partial sum of "pos" is achieved, where weight-0 nodes (e.g. dummy nodes) are skipped.
     */
    void searchPos
    (
     uint64_t & pos, //!< [in,out] Give global position to search. It is modified to relative position in bottom node.
     BTreeNodeT *& retNode, //!< [out] To capture parent of returned bottom node.
     uint8_t & retIdx //!< [out] To capture sibling idx of returned bottom node.
     ) noexcept {
      assert(pos < this->getSumOfWeight());

      retNode = this;
      while (true) {
        retIdx = 0;
        auto array = retNode->getConstPtr_psum();
        while (pos >= array[retIdx + 1]) {
          ++retIdx;
        }
        pos -= array[retIdx];
        if (retNode->isBorder()) {
          return;
        }
        retNode = retNode->getChildPtr(retIdx);
      }
    }


  private:
    //////////////////////////////// private modifier (intend to call them from member function of BTreeNode)
    void overflowToL
    (
     BTreeNodeT * lnode,
     BTreeNodeT * sndHalf,
     const uint8_t idx
     ) {
      assert(idx < numChildren_);

      const uint8_t lnum = lnode->getNumChildren();
      const uint8_t numToLeft = lnum/2 + kB/2 - lnum;
      const bool idxIsInLeft = idx < numToLeft;
      { // update left node (possibly first part)
        const uint8_t stop = (idxIsInLeft)? idx + 1 : numToLeft;
        for (uint8_t i = 0; i < stop; ++i) {
          lnode->pushbackBTreeNode(children_[i]);
        }
      }
      if (idxIsInLeft) { // update left node (second part)
        lnode->pushbackBTreeNode(sndHalf);
        for (uint8_t i = idx + 1; i < numToLeft; ++i) {
          lnode->pushbackBTreeNode(children_[i]);
        }
      }

      { // update this node (possibly first part)
        const uint8_t stop = (idxIsInLeft)? numChildren_ : idx+1;
        for (uint8_t i = numToLeft; i < stop; ++i) {
          putBTreeNode(children_[i], i - numToLeft);
        }
      }
      if (!idxIsInLeft) { // update this node (second part)
        putBTreeNode(sndHalf, idx + 1 - numToLeft);
        for (uint8_t i = idx + 1; i < numChildren_; ++i) {
          putBTreeNode(children_[i], i - numToLeft + 1);
        }
      }
      numChildren_ -= (numToLeft - !idxIsInLeft);

      //// update lmJumpNode
      if (!(this->isBorder())) {
        this->setLmJumpNode(this->children_[0]->getLmJumpNode());
      } else if (this->isJumpToBtm()) {
        this->setLmJumpNode(this->children_[0]);
      }
    }


    void overflowToR
    (
     BTreeNodeT * rnode,
     BTreeNodeT * sndHalf,
     const uint8_t idx
     ) {
      assert(numChildren_ == kB);
      assert(idx < kB);

      const uint8_t rnum = rnode->getNumChildren();
      const uint8_t leftEnd = rnum/2 + kB/2;
      const bool idxIsInLeft = idx < leftEnd;
      {
        const uint8_t shift = kB - leftEnd + !idxIsInLeft;
        rnode->shiftR(psum_[kB] - psum_[leftEnd], 0, shift); // Shift elements in rnode
      }

      numChildren_ = leftEnd;
      if (idxIsInLeft) {
        for (uint8_t i = leftEnd; i < kB; ++i) {
          rnode->putBTreeNode(children_[i], i - leftEnd);
        }
        psum_[idx + 1] -= sndHalf->getSumOfWeight();
        shiftR(0, idx + 1, 1);
        putBTreeNode(sndHalf, idx + 1);
      } else {
        for (uint8_t i = leftEnd; i < idx + 1; ++i) {
          rnode->putBTreeNode(children_[i], i - leftEnd);
        }
        rnode->putBTreeNode(sndHalf, idx + 1 - leftEnd);
        for (uint8_t i = idx + 1; i < kB; ++i) {
          rnode->putBTreeNode(children_[i], i + 1 - leftEnd);
        }
      }

      //// update lmJumpNode
      if (!(rnode->isBorder())) {
        rnode->setLmJumpNode(rnode->children_[0]->getLmJumpNode());
      } else if (rnode->isJumpToBtm()) {
        rnode->setLmJumpNode(rnode->children_[0]);
      }
    }


    void overflowToL_Btm
    (
     BTreeNodeT * lnode,
     BTreeNodeT * sndHalf,
     const uint64_t weight,
     const uint8_t idx
     ) {
      assert(idx < numChildren_);

      uint64_t temp_weight;
      const uint8_t lnum = lnode->getNumChildren();
      const uint8_t numToLeft = lnum/2 + kB/2 - lnum;
      const bool idxIsInLeft = idx < numToLeft;
      { // update lnode (possibly first part)
        const auto stop = (idxIsInLeft)? idx : numToLeft;
        for (uint8_t i = 0; i < stop; ++i) {
          temp_weight = getWeightOfChild(i);
          lnode->pushbackBtm(children_[i], temp_weight);
        }
      }
      if (idxIsInLeft) { // update lnode (second part)
        temp_weight = getWeightOfChild(idx) - weight;
        lnode->pushbackBtm(children_[idx], temp_weight);
        lnode->pushbackBtm(sndHalf, weight);
        for (uint8_t i = idx + 1; i < numToLeft; ++i) {
          temp_weight = getWeightOfChild(i);
          lnode->pushbackBtm(children_[i], temp_weight);
        }
      }

      { // update this node (possibly first part)
        const auto stop = (idxIsInLeft)? numChildren_ : idx;
        for (uint8_t i = numToLeft; i < stop; ++i) {
          temp_weight = getWeightOfChild(i);
          putBtm(children_[i], i - numToLeft, temp_weight);
        }
      }
      if (!idxIsInLeft) { // update this node (second part)
        temp_weight = getWeightOfChild(idx) - weight;
        putBtm(children_[idx], idx - numToLeft, temp_weight);
        if (numToLeft != 1) {
          putBtm(sndHalf, idx + 1 - numToLeft, weight);
          for (uint8_t i = idx + 1; i < numChildren_; ++i) {
            temp_weight = getWeightOfChild(i);
            putBtm(children_[i], i - numToLeft + 1, temp_weight);
          }
        } else { // need special treatment when numToLeft == 1
          temp_weight = psum_[idx + 1] - (psum_[idx] + weight); // This amount is equivalent to original psum_[1] (that is moved to lnode)
          putBtm(sndHalf, idx, weight);
          for (uint8_t i = idx + 2; i <= numChildren_; ++i) {
            psum_[i] -= temp_weight;
          }
        }
      }
      numChildren_ -= (numToLeft - !idxIsInLeft);

      //// update lmJumpNode
      if (this->isJumpToBtm()) {
        this->setLmJumpNode(this->children_[0]);
      }
    }


    void overflowToR_Btm
    (
     BTreeNodeT * rnode,
     BTreeNodeT * sndHalf,
     const uint64_t weight,
     const uint8_t idx
     ) {
      assert(numChildren_ == kB);
      assert(idx < kB);

      uint64_t temp_weight;

      const uint8_t rnum = rnode->getNumChildren();
      const uint8_t leftEnd = rnum/2 + kB/2;
      const bool idxIsInLeft = idx < leftEnd;
      { // shift elements in rnode
        const uint8_t shift = kB - leftEnd + !idxIsInLeft;
        temp_weight = psum_[kB] - psum_[leftEnd];
        rnode->shiftR_Btm(temp_weight, 0, shift);
      }
      { // move elements in this node to rnode (possibly first part)
        const auto stop = (idxIsInLeft)? kB : idx;
        for (uint8_t i = leftEnd; i < stop; ++i) {
          temp_weight = getWeightOfChild(i);
          rnode->putBtm(children_[i], i - leftEnd, temp_weight);
        }
      }
      if (!idxIsInLeft) { // move elements in this node to rnode (second part)
        temp_weight = getWeightOfChild(idx) - weight;
        rnode->putBtm(children_[idx], idx - leftEnd, temp_weight);
        rnode->putBtm(sndHalf, idx + 1 - leftEnd, weight);
        for (uint8_t i = idx + 1; i < kB; ++i) {
          temp_weight = getWeightOfChild(i);
          rnode->putBtm(children_[i], i + 1 - leftEnd, temp_weight);
        }
      }

      { // update this node
        numChildren_ = leftEnd;
        if (idxIsInLeft) {
          psum_[idx + 1] -= weight;
          shiftR_Btm(0, idx + 1, 1);
          putBtm(sndHalf, idx + 1, weight);
        }
      }

      //// update lmJumpNode
      if (rnode->isJumpToBtm()) {
        rnode->setLmJumpNode(rnode->children_[0]);
      }
    }


  public:
    //////////////////////////////// public modifier
    void setLmJumpNode
    (
     void * lmJumpNode
     ) noexcept {
      lmJumpNode_ = reinterpret_cast<BTreeNodeT *>(lmJumpNode);
    }


    void updateLmJumpNode
    (
     void * lmJumpNode
     ) noexcept {
      auto * node = this;
      while (true) {
        node->setLmJumpNode(lmJumpNode);
        if (node->isRoot() || node->getIdxInSibling() > 0) {
          break;
        }
        node = node->getParent();
      }
    }


    void unroot() noexcept {
      flags_ &= ~isRootBit;
    }


    void setParentRef
    (
     BTreeNodeT * newParent,
     uint8_t newIdxInSibling
     ) noexcept {
      this->parent_ = newParent;
      this->idxInSibling_ = newIdxInSibling;
    }


    void setChildPtr
    (
     BTreeNodeT * child,
     uint8_t idx
     ) noexcept {
      assert(idx < numChildren_);

      children_[idx] = child;
    }


    void makeNewRoot
    (
     BTreeNodeT * sndHalf
     ) noexcept {
      auto newRoot = new BTreeNodeT(this->getLmJumpNode(), true, false, this->isJumpToBtm(), this->isUnderSuperRoot());
      auto * parent = this->getParent();
      if (parent != nullptr) {
        if (this->isUnderSuperRoot()) { // parent is super root
          reinterpret_cast<SuperRootT *>(parent)->setRoot(newRoot);
        } else { // BTrees are stacked
          const auto idxInSib = this->getIdxInSibling();
          parent->setChildPtr(newRoot, idxInSib); // parent points to newRoot
          newRoot->setParentRef(parent, idxInSib); // newRoot points to parent
        }
      }
      newRoot->pushbackBTreeNode(this);
      newRoot->pushbackBTreeNode(sndHalf);
    }


    /*!
     * @brief Put node of ::BTreeNodeT type at idx.
     * @note
     *   - "psum_" is set using immediately preceeding psum value and weights of "child".
     *   - Links between "this" and "child" is fixed.
     *   - "numChildren_" is not incremented.
     */
    void putBTreeNode
    (
     BTreeNodeT * child,
     const uint8_t idx
     ) noexcept {
      assert(idx < kB);

      children_[idx] = child;
      psum_[idx + 1] = psum_[idx] + child->getSumOfWeight();
      child->setParentRef(this, idx);
    }


    /*!
     * @brief Pushback node of ::BTreeNodeT type.
     * @note
     *   - Links between "this" and "child" is fixed.
     *   - "numChildren_" is incremented.
     */
    void pushbackBTreeNode
    (
     BTreeNodeT * child
     ) noexcept {
      assert(numChildren_ < kB);

      putBTreeNode(child, numChildren_);
      ++numChildren_;
    }


    /*!
     * @brief Put bottom node other than ::BTreeNodeT type.
     * @note
     *   Link from "child" to "this" should be set outside this function (if "child" maintains it).
     */
    void putBtm
    (
     void * child,
     const uint8_t idx,
     const uint64_t weight
     ) noexcept {
      assert(isBorder());
      assert(idx < kB);

      children_[idx] = reinterpret_cast<BTreeNodeT *>(child);
      psum_[idx + 1] = psum_[idx] + weight;
    }


    /*!
     * @brief Pushback bottom node other than ::BTreeNodeT type.
     * @note
     *   Link from "child" to "this" should be set outside this function (if "child" maintains it).
     */
    void pushbackBtm
    (
     void * child,
     const uint64_t weight
     ) noexcept {
      assert(isBorder());
      assert(numChildren_ < kB);

      putBtm(child, numChildren_, weight);
      ++numChildren_;
    }


    void shiftR
    (
     const uint64_t add_weight,
     const uint8_t idx, //!< Beginning idx to shift.
     const uint8_t shift
     ) noexcept {
      assert(idx <= numChildren_);
      assert(numChildren_ + shift <= kB);

      for (uint8_t i = numChildren_; idx < i; --i) {
        auto child = children_[i-1];
        children_[i + shift - 1] = child;
        child->idxInSibling_ = i + shift - 1;
      }
      for (uint8_t i = numChildren_; idx < i; --i) {
        psum_[i + shift] = psum_[i] + add_weight;
      }
      numChildren_ += shift;
    }


    void shiftR_Btm
    (
     const uint64_t add_weight,
     const uint8_t idx, //!< Beginning idx to shift.
     const uint8_t shift
     ) noexcept {
      assert(idx <= numChildren_);
      assert(numChildren_ + shift <= kB);

      for (uint8_t i = numChildren_; idx < i; --i) {
        children_[i + shift - 1] = children_[i-1];
      }
      for (uint8_t i = numChildren_; idx < i; --i) {
        psum_[i + shift] = psum_[i] + add_weight;
      }
      numChildren_ += shift;
    }


    /*!
     * @brief Handle the situation where 'children_[idx]' is split to 'children_[idx]' and 'sndHalf' when child node is ::BTreeNode type.
     */
    void handleSplitOfChild
    (
     BTreeNodeT * sndHalf,
     const uint8_t idx
     ) {
      assert(idx <= numChildren_);

      const auto end = numChildren_;
      if (end < kB) { // Easy case: Current node is not full.
        numChildren_ = idx;
        this->pushbackBTreeNode(children_[idx]);
        auto * pushC = sndHalf;
        for (uint8_t i = idx+1; i <= end; ++i) {
          auto * tmp = children_[i];
          this->pushbackBTreeNode(pushC);
          pushC = tmp;
        }
        return;
      }

      if (!isRoot()) {
        //// Check siblings if they have space.
        if (idxInSibling_) { // Check previous sibling.
          auto lnode = parent_->getChildPtr(idxInSibling_ - 1);
          if (lnode->getNumChildren() < kB - 1) { // Previous sibling is not full. (-1 for easier implementation)
            overflowToL(lnode, sndHalf, idx);
            if (!(this->isRoot())) {
              const uint64_t psum = parent_->getPSum(idxInSibling_ - 1) + lnode->getSumOfWeight(); // psum value for lnode in its parent
              parent_->changePSumAt(idxInSibling_ - 1, psum);
            }
            return;
          }
        }
        if (idxInSibling_ + 1 < parent_->getNumChildren()) { // Check next sibling.
          auto rnode = parent_->getChildPtr(idxInSibling_ + 1);
          if (rnode->getNumChildren() < kB - 1) { // Next sibling is not full. (-1 for easier implementation)
            overflowToR(rnode, sndHalf, idx);
            if (!(this->isRoot())) {
              const uint64_t psum = parent_->getPSum(idxInSibling_) + this->getSumOfWeight(); // psum value for this node in its parent
              parent_->changePSumAt(idxInSibling_, psum);
            }
            return;
          }
        }
      }

      { // this node has to be split
        auto newNode = new BTreeNodeT(nullptr, false, this->isBorder(), this->isJumpToBtm(), this->isUnderSuperRoot());
        overflowToR(newNode, sndHalf, idx);
        if (!isRoot()) {
          parent_->handleSplitOfChild(newNode, idxInSibling_);
        } else {
          this->unroot();
          this->makeNewRoot(newNode);
        }
      }
    }



    /*!
     * @brief Merge children of rnode into this node
     * @attention children of "this" node should be BTreeNodeT type.
     */
    void mergeBTreeChildrenWithOneDelete
    (
     BTreeNodeT * rnode,
     const uint8_t idx,
     const bool isIdxOfLeft
     ) {
      assert(static_cast<uint16_t>(idx) < static_cast<uint16_t>(kB) * 2);

      const auto rnum = rnode->getNumChildren();
      if (isIdxOfLeft) {
        const auto lnum = this->getNumChildren();
        for (uint8_t i = idx + 1; i < lnum; ++i) {
          putBTreeNode(children_[i], i - 1);
        }
        --numChildren_;
        for (uint8_t i = 0; i < rnum; ++i) {
          pushbackBTreeNode(rnode->getChildPtr(i));
        }
      } else {
        for (uint8_t i = 0; i < idx; ++i) {
          pushbackBTreeNode(rnode->getChildPtr(i));
        }
        for (uint8_t i = idx + 1; i < rnum; ++i) {
          pushbackBTreeNode(rnode->getChildPtr(i));
        }
      }
    }


    /*!
     * @brief Handle deletion of child of btree node type at idx
     * @attention assume that tree will not be empty
     */
    void handleDeleteOfChild
    (
     const uint8_t idx
     ) {
      // {//debug
      //   std::cerr << __func__ << this << ": idx = " << (int)idx << std::endl;
      // }
      assert(idx < numChildren_);

      const uint8_t num_new = numChildren_ - 1;
      if (num_new == 0) {
        parent_->handleDeleteOfChild(idxInSibling_);
        delete this;
        return;
      }

      if (idx == 0) {
        this->updateLmJumpNode(this->getChildPtr(1)->getLmJumpNode());
      }

      if (!isRoot() && num_new < kB / 2) {
        if (idxInSibling_) { // Check previous sibling if it can be merged
          auto lnode = parent_->getChildPtr(idxInSibling_ - 1);
          const uint16_t lnum = lnode->getNumChildren();
          if (lnum + num_new <= 3 * kB / 4) {
            lnode->mergeBTreeChildrenWithOneDelete(this, idx, false);
            parent_->changePSumAt(idxInSibling_ - 1, parent_->getPSum(idxInSibling_ + 1));
            parent_->handleDeleteOfChild(idxInSibling_);
            delete this;
            return;
          }
        }
        if (idxInSibling_ + 1 < parent_->getNumChildren()) { // Check next sibling if it can be merged
          auto rnode = parent_->getChildPtr(idxInSibling_ + 1);
          const uint16_t rnum = rnode->getNumChildren();
          if (rnum + num_new <= 3 * kB / 4) {
            this->mergeBTreeChildrenWithOneDelete(rnode, idx, true);
            parent_->changePSumAt(idxInSibling_, parent_->getPSum(idxInSibling_ + 2));
            parent_->handleDeleteOfChild(idxInSibling_ + 1);
            delete rnode;
            return;
          }
        }
      }

      // {//debug
      //   std::cout << __FUNCTION__ << ": easy case" << std::endl;
      // }
      for (uint8_t i = idx + 1; i < numChildren_; ++i) {
        putBTreeNode(children_[i], i - 1);
      }
      numChildren_ = num_new;

      if (isRoot() && num_new == 1 && parent_ != nullptr) {
        BTreeNodeT * newRoot = children_[0];
        newRoot->setParentRef(parent_, idxInSibling_);
        newRoot->flags_ |= isRootBit; // enRoot
        if (this->isUnderSuperRoot()) {
          reinterpret_cast<SuperRootT *>(parent_)->setRoot(newRoot);
        } else {
          parent_->setChildPtr(newRoot, idxInSibling_);
        }
        delete this;
      }
    }


    /*!
     * @brief Handle the situation where 'children_[idx]' is split to 'children_[idx]' and 'sndHalf' when child node is not ::BTreeNode type.
     */
    uint8_t handleSplitOfBtm
    (
     BTreeNodeT * sndHalf,
     const uint64_t weight,
     const uint8_t idx
     ) {
      assert(isBorder());
      assert(idx <= numChildren_);

      if (numChildren_ < kB) { // Easy case: Current node is not full.
        psum_[idx + 1] -= weight;
        shiftR_Btm(0, idx + 1, 1);
        putBtm(sndHalf, idx + 1, weight);
        return 0;
      }

      if (!isRoot()) {
        // Check siblings if they have space.
        if (idxInSibling_) { // Check previous sibling.
          auto lnode = parent_->getChildPtr(idxInSibling_ - 1);
          const uint8_t lnum = lnode->getNumChildren();
          if (lnum < kB - 1) { // Previous sibling is not full. (-1 for easier implementation)
            overflowToL_Btm(lnode, sndHalf, weight, idx);
            if (!(this->isRoot())) {
              const uint64_t psum = parent_->getPSum(idxInSibling_ - 1) + lnode->getSumOfWeight(); // psum value for lnode in its parent
              parent_->changePSumAt(idxInSibling_ - 1, psum);
            }
            return lnum/2 + kB/2 - lnum; // return "numToLeft" in overflowTL_Btm. Actual number of increase of lnode is "numToLeft + (idx < numToLeft)".
          }
        }
        if (idxInSibling_ + 1 < parent_->getNumChildren()) { // Check next sibling.
          auto rnode = parent_->getChildPtr(idxInSibling_ + 1);
          if (rnode->getNumChildren() < kB - 1) { // Next sibling is not full. (-1 for easier implementation)
            overflowToR_Btm(rnode, sndHalf, weight, idx);
            if (!(this->isRoot())) {
              const uint64_t psum = parent_->getPSum(idxInSibling_) + this->getSumOfWeight(); // psum value for this node in its parent
              parent_->changePSumAt(idxInSibling_, psum);
            }
            return 0;
          }
        }
      }

      { // This node has to be split.
        auto newNode = new BTreeNode(children_[kB/2], false, true, this->isJumpToBtm(), this->isUnderSuperRoot());
        overflowToR_Btm(newNode, sndHalf, weight, idx);
        if (!isRoot()) {
          parent_->handleSplitOfChild(newNode, idxInSibling_);
        } else {
          this->unroot();
          this->makeNewRoot(newNode);
        }
        return 0;
      }
    }


    /*!
     * @brief Merge children of rnode into this node
     * @attention children of "this" node should be BTreeNodeT type.
     */
    void mergeBtmChildrenWithOneDelete
    (
     BTreeNodeT * rnode,
     const uint8_t idx,
     const bool isIdxOfLeft
     ) {
      // {//debug
      //   std::cerr << __func__ << "(" << this << "): rnode = " << rnode
      //             << ", idx = " << (int)idx << ", isIdxOfLeft = " << isIdxOfLeft << std::endl;
      // }
      assert(isBorder());
      assert(static_cast<uint16_t>(idx) < static_cast<uint16_t>(kB) * 2);

      const auto rnum = rnode->getNumChildren();
      if (isIdxOfLeft) {
        const auto lnum = this->getNumChildren();
        for (uint8_t i = idx + 1; i < lnum; ++i) {
          putBtm(children_[i], i - 1, getWeightOfChild(i));
        }
        --numChildren_;
        for (uint8_t i = 0; i < rnum; ++i) {
          pushbackBtm(rnode->getChildPtr(i), rnode->getWeightOfChild(i));
        }
      } else {
        for (uint8_t i = 0; i < idx; ++i) {
          pushbackBtm(rnode->getChildPtr(i), rnode->getWeightOfChild(i));
        }
        for (uint8_t i = idx + 1; i < rnum; ++i) {
          pushbackBtm(rnode->getChildPtr(i), rnode->getWeightOfChild(i));
        }
      }
    }


    /*!
     * @brief Handle deletion of bottom child at idx
     */
    // void handleDeleteOfBtm
    // (
    //  uint8_t & idx, //!< [in, out]
    //  BTreeNodeT * retParent //!< [out]
    //  ) {
    //   {//debug
    //     std::cerr << __func__ << this << ": idx = " << (int)idx << ", numChildren_ = " << (int)numChildren_ << std::endl;
    //   }
    //   assert(isBorder());
    //   assert(idx < numChildren_);

    //   retParent = this;
    //   const uint8_t num_new = numChildren_ - 1;
    //   if (num_new == 0) {
    //     retParent = getPrevBtmRef(idx); // idx is modified
    //     parent_->handleDeleteOfChild(idxInSibling_);
    //     delete this;
    //     return;
    //   }
    //   if (!isRoot() && num_new < kB / 2) {
    //     if (idxInSibling_) { // Check previous sibling if it can be merged
    //       auto lnode = parent_->getChildPtr(idxInSibling_ - 1);
    //       const uint16_t lnum = lnode->getNumChildren();
    //       if (lnum + num_new <= 3 * kB / 4) {
    //         lnode->mergeBtmChildrenWithOneDelete(this, idx, false);
    //         parent_->changePSumAt(idxInSibling_ - 1, parent_->getPSum(idxInSibling_ + 1));
    //         parent_->handleDeleteOfChild(idxInSibling_);
    //         delete this;
    //         idx = static_cast<uint8_t>(lnum + idx - 1);
    //         retParent = lnode;
    //         return;
    //       }
    //     }
    //     if (idxInSibling_ + 1 < parent_->getNumChildren()) { // Check next sibling if it can be merged
    //       auto rnode = parent_->getChildPtr(idxInSibling_ + 1);
    //       const uint16_t rnum = rnode->getNumChildren();
    //       if (rnum + num_new <= 3 * kB / 4) {
    //         this->mergeBtmChildrenWithOneDelete(rnode, idx, true);
    //         parent_->changePSumAt(idxInSibling_, parent_->getPSum(idxInSibling_ + 2));
    //         parent_->handleDeleteOfChild(idxInSibling_ + 1);
    //         delete rnode;
    //         --idx;
    //         return;
    //       }
    //     }
    //   }

    //   // {//debug
    //   //   std::cout << __FUNCTION__ << ": easy case" << std::endl;
    //   // }
    //   for (uint8_t i = idx + 1; i < numChildren_; ++i) {
    //     putBtm(children_[i], i - 1, getWeightOfChild(i));
    //   }
    //   numChildren_ = num_new;

    //   if (isRoot() && num_new == 1 && parent_ != nullptr) {
    //     BTreeNodeT * newRoot = children_[0];
    //     newRoot->setParentRef(parent_, idxInSibling_);
    //     if (this->isUnderSuperRoot()) {
    //       reinterpret_cast<SuperRootT *>(parent_)->setRoot(newRoot);
    //     } else {
    //       parent_->setChildPtr(newRoot, idxInSibling_);
    //     }
    //     idx = 0;
    //     retParent = newRoot;
    //   } else {
    //     --idx;
    //   }
    // }


    /*!
     * @brief Handle deletion of bottom child at idx
     */
    void handleDeleteOfBtm
    (
     uint8_t idx //!<
     ) {
      // {//debug
      //   std::cerr << __func__ << "(" << this << "): idx = " << (int)idx << ", numChildren_ = " << (int)numChildren_ << std::endl;
      // }
      assert(isBorder());
      assert(idx < numChildren_);

      const uint8_t num_new = numChildren_ - 1;
      if (num_new == 0) {
        parent_->handleDeleteOfChild(idxInSibling_);
        delete this;
        return;
      }

      if (idx == 0 && this->isJumpToBtm()) {
        this->updateLmJumpNode(this->getChildPtr(1));
      }

      if (!isRoot() && num_new < kB / 2) {
        if (idxInSibling_) { // Check previous sibling if it can be merged
          auto lnode = parent_->getChildPtr(idxInSibling_ - 1);
          const uint16_t lnum = lnode->getNumChildren();
          if (lnum + num_new <= 3 * kB / 4) {
            lnode->mergeBtmChildrenWithOneDelete(this, idx, false);
            parent_->changePSumAt(idxInSibling_ - 1, parent_->getPSum(idxInSibling_ + 1));
            parent_->handleDeleteOfChild(idxInSibling_);
            delete this;
            return;
          }
        }
        if (idxInSibling_ + 1 < parent_->getNumChildren()) { // Check next sibling if it can be merged
          auto rnode = parent_->getChildPtr(idxInSibling_ + 1);
          const uint16_t rnum = rnode->getNumChildren();
          if (rnum + num_new <= 3 * kB / 4) {
            this->mergeBtmChildrenWithOneDelete(rnode, idx, true);
            parent_->changePSumAt(idxInSibling_, parent_->getPSum(idxInSibling_ + 2));
            parent_->handleDeleteOfChild(idxInSibling_ + 1);
            delete rnode;
            return;
          }
        }
      }

      // {//debug
      //   std::cerr << __func__ << "(" << this << "): easy case: idx = " << (int)idx << std::endl;
      //   if (idxInSibling_ == 6 && numChildren_ == 11 && idx == 0) {
      //     this->printStatistics(std::cerr, false);
      //   }
      // }
      for (uint8_t i = idx + 1; i < numChildren_; ++i) {
        putBtm(children_[i], i - 1, getWeightOfChild(i));
      }
      numChildren_ = num_new;
    }


    /*!
     * @brief Change weight of "idx"-th child of this node, and accordingly all "psum_" values needed to be fixed.
     */
    void changePSumFrom
    (
     const uint8_t idx,
     const int64_t change
     ) noexcept {
      assert(idx < numChildren_);

      for (uint8_t i = idx; i < numChildren_; ++i) {
        psum_[i+1] += static_cast<uint64_t>(change);
      }
      if (!(isRoot() && isUnderSuperRoot()) && parent_ != nullptr) { // When we stack two or more BTrees, the change will be propagated.
        parent_->changePSumFrom(idxInSibling_, change);
      }
    }


    /*!
     * @brief Change psum value for "idx"-th child.
     */
    void changePSumAt
    (
     const uint8_t idx,
     const uint64_t psum_value
     ) noexcept {
      assert(idx < numChildren_);

      psum_[idx + 1] = psum_value;
    }


    /*!
     * @brief Pushback bottom node other than ::BTreeNodeT type.
     * @note
     *   Link from "child" to "this" should be set outside this function (if "child" maintains it).
     */
    void putFirstBtm
    (
     void * child,
     const uint64_t weight
     ) noexcept {
      assert(isBorder());
      assert(numChildren_ == 0);

      ++numChildren_;
      putBtm(child, 0, weight);
      if (isJumpToBtm()) {
        setLmJumpNode(child);
      }
    }


  public:
    //// calculate statistics
    size_t calcMemBytes
    (
     bool includeThis = true
     ) const noexcept {
      size_t sumOfBytes = sizeof(*this) * includeThis;
      if (!isBorder()) {
        for (uint8_t i = 0; i < numChildren_; ++i) {
          sumOfBytes += children_[i]->calcMemBytes();
        }
      }
      return sumOfBytes;
    }


    size_t calcNumUsed() const noexcept {
      size_t numOfUsed = numChildren_;
      if (!isBorder()) {
        for (uint8_t i = 0; i < numChildren_; ++i) {
          numOfUsed += children_[i]->calcNumUsed();
        }
      }
      return numOfUsed;
    }


    size_t calcNumSlots() const noexcept {
      size_t numOfSlots = kB;
      if (!isBorder()) {
        for (uint8_t i = 0; i < numChildren_; ++i) {
          numOfSlots += children_[i]->calcNumSlots();
        }
      }
      return numOfSlots;
    }


    void printStatistics
    (
     std::ostream & os,
     const bool verbose = false
     ) const noexcept {
      os << "BTreeNode object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
      os << "parent = " << parent_ << ", idx in sibling = " << static_cast<uint64_t>(idxInSibling_)
         << ", num of children = " << static_cast<uint64_t>(numChildren_) << ", BTree arity param kB = " << static_cast<uint64_t>(kB) << std::endl;
      os << "leftmost jump node = " << lmJumpNode_
         << ", flags (isRoot, isBorder, isJumpToBtm, isUnderSuperRoot, isDummy) = ("
         << isRoot() << ", " << isBorder() << ", " << isJumpToBtm() << ", " << isUnderSuperRoot() << ", " << isDummy() << ")" << std::endl;
      os << "dump partial sums (child pointer): " << std::endl;
      os << getPSum(0) << ", ";
      for (uint8_t i = 0; i < numChildren_; ++i) {
        os << "[" << static_cast<uint64_t>(i) << "] " << getPSum(i+1) << " (" << getChildPtr(i) << "), ";
      }
      os << std::endl;
      os << "BTreeNode object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
      if (verbose && !isBorder()) {
        for (uint8_t i = 0; i < numChildren_; ++i) {
          getChildPtr(i)->printStatistics(os, verbose);
        }
      }
    }
  };
} // namespace itmmti

#endif
