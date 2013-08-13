      integer function blocklow(id,p,n)
      integer id,p,n
      blocklow = id * n / p
      return
      end

      integer function blockhigh(id,p,n)
      integer id,p,n
      integer blocklow
      blockhigh = blocklow(id + 1, p, n) - 1
      return
      end

      integer function blocksize(id,p,n)
      integer id,p,n
      integer blocklow, blockhigh
      blocksize = blocklow(id + 1,p,n) - blocklow(id,p,n)
      return
      end

      integer function blockowner(bindex,p,n)
      integer bindex,p,n
      blockowner = (p * (bindex + 1) - 1) / n
      return
      end
